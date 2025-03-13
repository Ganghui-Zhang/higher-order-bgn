#pragma once
#include <cmath>
#include <type_traits>
#include <vector>

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/common/parametertree.hh>

#include <dune/curvedgrid/curvedgrid.hh>
#include <dune/curvedgrid/gridfunctions/discretegridviewfunction.hh>

#include <dune/functions/functionspacebases/defaultglobalbasis.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/functions/functionspacebases/compositebasis.hh>
#include <dune/functions/functionspacebases/subspacebasis.hh>
#include <dune/functions/functionspacebases/interpolate.hh>
#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>

#include <dune/geometry/quadraturerules.hh>

#include <dune/istl/bvector.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/matrix.hh>
#include <dune/istl/matrixindexset.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/umfpack.hh>

#include <dune/vtk/pvdwriter.hh>
#include <dune/vtk/vtkwriter.hh>
#include <dune/vtk/datacollectors/lagrangedatacollector.hh>

class LocalAssembler
{
public:
  LocalAssembler(Dune::ParameterTree const& pt)
    : tau_(pt.get<double>("adapt.timestep",0.1))
    , quadOrder_(pt.get<int>("solution.quad_order",4))
  {}

  // Compute the stiffness matrix for a single element
  template <class LocalView, class LocalParametrization,
            class Geometry, class MatrixType, class VectorType>
  void operator()(LocalView const& localView, LocalParametrization const& X_e,
                  Geometry const& geometry, MatrixType& elementMatrix, VectorType& elementVector) const
  {
    auto const& element = localView.element();
    auto const& node = localView.tree();

    // Total number of DOFs on the element
    std::size_t s = node.size();

    using namespace Dune::Indices;

    auto const& PositionlocalFE = node.child(_0,0).finiteElement();
    auto const& CurvaturelocalFE = node.child(_1).finiteElement();
    auto const& localBasis = CurvaturelocalFE.localBasis();

    // Number of DOFs on the element per component
    std::size_t s0 = localBasis.size();

    // Set all matrix entries to zero
    elementMatrix.setSize(s,s);
    elementMatrix = 0;

    // Set all vector entries to zero
    elementVector.resize(s);
    elementVector = 0;

    using GlobalCoordinate = typename Geometry::GlobalCoordinate;
    using Traits = typename std::decay_t<decltype(localBasis)>::Traits;
    std::vector<typename Traits::RangeType> shapeValues;
    std::vector<typename Traits::JacobianType> shapeGradients;
    std::vector<GlobalCoordinate> gradients;

    // Get a quadrature rule
    const auto& quad = Dune::QuadratureRules<double,Geometry::mydimension>::rule(element.type(), quadOrder_);

    // Loop over all quadrature points
    for (const auto& qp : quad)
    {
      // integration element dG_{s-1}
      const auto dS = geometry.integrationElement(qp.position()) * qp.weight();

      // the inverse of the transposed geometry Jacobian
      const auto& Jtinv = geometry.jacobianInverseTransposed(qp.position());

      auto const& n = geometry.normal(qp.position());

      // Evaluate all shape function values, normal and gradients at quadrature point
      localBasis.evaluateFunction(qp.position(), shapeValues);
      localBasis.evaluateJacobian(qp.position(), shapeGradients);

      // Compute the shape function gradients on the real element
      gradients.resize(shapeGradients.size());
      for (std::size_t j = 0; j < gradients.size(); ++j)
        Jtinv.mv(shapeGradients[j][0], gradients[j]);


      // Compute the actual matrix entries
      for (std::size_t j1 = 0; j1 < s0; ++j1) {
        for (std::size_t j2 = 0; j2 < s0; ++j2) {
          const auto mass = shapeValues[j1] * shapeValues[j2] * dS;
          const auto laplace = dot(gradients[j1], gradients[j2]) * dS;

          std::size_t row = node.child(_1).localIndex(j1);
          std::size_t col = node.child(_1).localIndex(j2);

          // tau * (v,w) for mcf
           elementMatrix[row][col] += - tau_ * mass;

          // tau * (grad(v), grad(w)) for sd
          // elementMatrix[row][col] += -tau_ * laplace;

          for (int i = 0; i < n.size(); ++i) {
            const auto mass_normal = shapeValues[j1] * shapeValues[j2] * n[i] * dS;
            std::size_t row_i = node.child(_0,i).localIndex(j1);
            std::size_t col_i = node.child(_0,i).localIndex(j2);

            // (v, normal*w)
            elementMatrix[row][col_i] += mass_normal;
            elementMatrix[col_i][row] += mass_normal;

            // (grad(v), grad(w))
            elementMatrix[row_i][col_i] += laplace;
          } // end i
        } // end j2
      } // end j1

      // Compute the actual vector entries
      Dune::FieldVector<double,2> fAtQP = X_e(qp.position());

      for (std::size_t j1 = 0; j1 < s0; ++j1) {
        for (int i = 0; i < n.size(); ++i) {
          std::size_t row = node.child(_1).localIndex(j1);

          // (x, normal*w)
          elementVector[row] += fAtQP[i] * shapeValues[j1] * n[i] * dS;
        } // end i
      } // end j1
    }
  }

private:
  double tau_;
  int quadOrder_;
};


struct Assembler
{
  // assembler matrix and vector on the grid
  template <class Basis, class Parametrization, class LocalAssembler, class BM, class BV>
  void operator()(Basis const& basis, Parametrization const& X, LocalAssembler const& localAssembler,
                  Dune::BCRSMatrix<BM>& matrix, Dune::BlockVector<BV>& rhs) const
  {
    // resize matrix and vector and initialize nonzero structure
    this->initialize(basis, matrix);
    this->initialize(basis, rhs);

    Dune::Matrix<double> elementMatrix;
    Dune::BlockVector<double> elementVector;

    auto X_e = localFunction(X);
    auto localView = basis.localView();
    for(const auto& e : elements(basis.gridView()))
    {
      X_e.bind(e);
      localView.bind(e);

      // construct the curved geometry from the parametrization
      Dune::LocalFunctionGeometry curvedGeometry{referenceElement(e), X_e};

      localAssembler(localView, X_e, curvedGeometry, elementMatrix, elementVector);
      this->scatter(localView, elementMatrix, elementVector, matrix, rhs);
    }
  }

private:
  template <class FEBasis, class B>
  void initialize(FEBasis const& feBasis, Dune::BCRSMatrix<B>& matrix) const
  {
    auto n = feBasis.size();
    Dune::MatrixIndexSet nb(n,n);

    auto localView = feBasis.localView();
    for(const auto& e : elements(feBasis.gridView()))
    {
      localView.bind(e);

      for (std::size_t i=0; i<localView.tree().size(); i++) {
        for (std::size_t j=0; j<localView.tree().size(); j++) {
          auto iIdx = localView.index(i);
          auto jIdx = localView.index(j);

          nb.add(iIdx[0], jIdx[0]);
        }
      }
    }

    nb.exportIdx(matrix);
    matrix = 0;
  }

  template <class FEBasis, class B>
  void initialize(FEBasis const& feBasis, Dune::BlockVector<B>& vector) const
  {
    vector.resize(feBasis.size());
    vector = 0;
  }

  template <class LocalView, class BM, class BV>
  void scatter(LocalView const& localView,
              Dune::Matrix<double> const& elMatrix, Dune::BlockVector<double> const& elVector,
              Dune::BCRSMatrix<BM>& matrix, Dune::BlockVector<BV>& rhs) const
  {
    for (std::size_t i = 0; i < elMatrix.N(); ++i) {
      auto row = localView.index(i);
      for (std::size_t j = 0; j < elMatrix.M(); ++j) {
        auto col = localView.index(j);
        matrix[row[0]][col[0]] += elMatrix[i][j];
      }

      rhs[row[0]] += elVector[i];
    }
  }
};

// Compute geometric quantities: surface area and volumn
  template <class Basis, class Parametrization>
  auto PerimeterAndArea(Dune::ParameterTree const& pt, Basis const& basis, Parametrization const& X)
  {
    double Perimeter = 0;
    double Area = 0;

    auto X_e = localFunction(X);
    auto localView = basis.localView();

    int ko = pt.get<int>("grid.ko", 1);

    for (auto const& e : elements(basis.gridView()))
    {
      X_e.bind(e);
      localView.bind(e);

      // construct the curved geometry from the parametrization
      Dune::LocalFunctionGeometry geometry{referenceElement(e), X_e};

      auto const& quadRule = Dune::QuadratureRules<double,1>::rule(e.type(), ko);
      for (auto const& qp : quadRule)
       {
         const auto dS = geometry.integrationElement(qp.position()) * qp.weight();
         auto const& n = geometry.normal(qp.position());

        //  surface_area = int_{gridview} dx
         Perimeter += dS;
          
         Dune::FieldVector<double,2> fAtQP = X_e(qp.position());
         for (int i = 0; i < 2; ++i) {
          //  volume = int_{gridview} 1/3*X_e*normal dx
           Area += 1.0/2 * fAtQP[i] * n[i] * dS;
         }
      }
    }

    return std::tuple{Perimeter,Area};
  }


  // Compute geometric quantities: surface area and volumn
  template <class Basis, class Parametrization>
  auto Meshratio(Dune::ParameterTree const& pt, Basis const& basis, Parametrization const& X)
  {
    double Meshratio = 0;
    
    std::vector<double> MeshratioVector;
    auto X_e = localFunction(X);
    auto localView = basis.localView();

    int ko = pt.get<int>("grid.ko", 1);

    for (auto const& e : elements(basis.gridView()))
    {
      
      X_e.bind(e);
      localView.bind(e);
      double Meshratio_e = 0;

      // construct the curved geometry from the parametrization
      Dune::LocalFunctionGeometry geometry{referenceElement(e), X_e};

      auto const& quadRule = Dune::QuadratureRules<double,1>::rule(e.type(), ko);
      for (auto const& qp : quadRule)
       {
         const auto dS = geometry.integrationElement(qp.position()) * qp.weight();

        //  surface_area = int_{gridview} dx
         Meshratio_e += dS;
      }
      MeshratioVector.push_back(Meshratio_e);
    }
    double max = *max_element(MeshratioVector.begin(), MeshratioVector.end());
    double min = *min_element(MeshratioVector.begin(), MeshratioVector.end());
    
    Meshratio  = max/min;

    return Meshratio;
  }

// Compute geometric quantities: surface area and volumn
  template <class Basis, class Parametrization>
  auto Meshsize(Dune::ParameterTree const& pt, Basis const& basis, Parametrization const& X)
  {
    double Meshsize = 0;
    
    std::vector<double> MeshsizeVector;
    auto X_e = localFunction(X);
    auto localView = basis.localView();

    int ko = pt.get<int>("grid.ko", 1);

    for (auto const& e : elements(basis.gridView()))
    {
      
      X_e.bind(e);
      localView.bind(e);
      double Meshsize_e = 0;

      // construct the curved geometry from the parametrization
      Dune::LocalFunctionGeometry geometry{referenceElement(e), X_e};

      auto const& quadRule = Dune::QuadratureRules<double,1>::rule(e.type(), ko);
      for (auto const& qp : quadRule)
       {
         const auto dS = geometry.integrationElement(qp.position()) * qp.weight();

        //  surface_area = int_{gridview} dx
         Meshsize_e += dS;
      }
      MeshsizeVector.push_back(Meshsize_e);
    }
    double max = *max_element(MeshsizeVector.begin(), MeshsizeVector.end());
    
    Meshsize  = max;

    return Meshsize;
  }

template <class HostGrid, class InitialSurface>
void run(Dune::ParameterTree const& pt, HostGrid& hostGrid, InitialSurface const& initialSurface)
{
  int kg = pt.get<int>("grid.kg", 1);
  int ko = pt.get<int>("grid.ko", 10);

  using namespace Dune::Functions::BasisFactory;
  auto feBasis = makeBasis(hostGrid.leafGridView(), 
    composite(
      power<HostGrid::dimensionworld>(lagrange(kg), flatInterleaved()),
      lagrange(kg),
      flatLexicographic()));

  auto positionBasis = Dune::Functions::subspaceBasis(feBasis, Dune::Indices::_0);
  auto curvatureBasis = Dune::Functions::subspaceBasis(feBasis, Dune::Indices::_1);

 
  using Matrix = Dune::BCRSMatrix<double>;
  Matrix matrix;
 
  using Vector = Dune::BlockVector<double>;
  Vector rhs;

  // interpolate initial surface to solution vector
  Vector solution;

  auto Xh = Dune::Functions::makeDiscreteGlobalBasisFunction<Dune::FieldVector<double,2>>(positionBasis, solution);
  auto Hh = Dune::Functions::makeDiscreteGlobalBasisFunction<double>(curvatureBasis, solution);

  Dune::Functions::interpolate(positionBasis, solution, initialSurface);

  // the computational domain
  // Dune::CurvedGrid grid{hostGrid, Xh};
  // auto gridView = grid.leafGridView();

  // // create a grid filewriter function to write solution of transformed grid
  // auto dataCollector = Dune::Vtk::LagrangeDataCollector{gridView, kg};
  // using Writer = decltype(Dune::Vtk::UnstructuredGridWriter{dataCollector});
  // auto pvdWriter = Dune::Vtk::PvdWriter<Writer>{dataCollector};
  // auto writeFiles = [&,filename=pt.get<std::string>("grid.output_filename", "grid")](double time)
  // {
  //   pvdWriter.writeTimestep(time, filename + ".vtu");
  // };

  Assembler assembler;
  LocalAssembler localAssembler{pt};

  double startTime = pt.get<double>("adapt.start_time", 0.0);
  double endTime = pt.get<double>("adapt.end_time", 0.03);
  double tau = pt.get<double>("adapt.timestep", 0.1);
  
  std::vector<double> time;
  // std::vector<double> perimeter;
  // std::vector<double> normalized_perimeter;
  // std::vector<double> area;
  // std::vector<double> relative_area_loss;
  // std::vector<double> meshratio;

  using std::sqrt;
  double tol = sqrt(std::numeric_limits<double>::epsilon());

  double t = startTime;
  int total_step = (endTime - startTime)/tau;
  //writeFiles(t);

  // auto[L0,A0] = PerimeterAndArea(pt,feBasis,Xh);
  // auto M0 = Meshratio(pt,feBasis,Xh);
  // auto h0 = Meshsize(pt,feBasis,Xh);


  // time.push_back(t);
  // perimeter.push_back(L0);
  // normalized_perimeter.push_back(L0/L0);
  // area.push_back(A0);
  // relative_area_loss.push_back((A0 - A0)/A0);
  // meshratio.push_back(M0);


  while (t + tau < endTime + tol)
  {
    t += tau;

    // assemble matrix and vector
    assembler(feBasis, Xh, localAssembler, matrix, rhs);

    // create linear solver
    auto solver = Dune::UMFPack<Matrix>(matrix, pt.get<int>("solver.verbose",0));

    // solve the linear system
    Dune::InverseOperatorResult statistics;
    solver.apply(solution, rhs, statistics);

    // write solution to file
    //writeFiles(t);
  
    // time.push_back(t);
    // perimeter.push_back(L);
    // normalized_perimeter.push_back(L/L0);
    // area.push_back(A);
    // relative_area_loss.push_back((A - A0)/A0);
    // meshratio.push_back(M);
  }
  
  // Output the Geometric quantities

  // std::ofstream time_file ("Time.txt");
  // if (time_file.is_open())
  // {
  //   for (int count = 0; count < total_step; count++)
  //       time_file  
  //       << std::fixed << std::setprecision(15) << time[count] << "\n";
  //   time_file.close();
  // }

  // std::ofstream Perimeter_file ("Normalized_Perimeter.txt");
  // if (Perimeter_file.is_open())
  // {
  //   for (int count = 0; count < total_step; count++)
  //       Perimeter_file  
  //       << std::fixed << std::setprecision(15) << normalized_perimeter[count] << "\n";
  //   // Close the file
  //   Perimeter_file.close();
  // }
    
  // std::ofstream Area_file ("Relative_area.txt");
  // if (Area_file.is_open())
  // {
  //   for (int count = 0; count < total_step; count++)
  //       Area_file  
  //       << std::fixed << std::setprecision(15) << relative_area_loss[count] << "\n";
  //   // Close the file
  //   Area_file.close();
  // }

  // std::ofstream Mesh_file ("Mesh_ratio.txt");
  // if (Mesh_file.is_open())
  // {
  //   for (int count = 0; count < total_step; count++)
  //       Mesh_file  
  //       << std::fixed << std::setprecision(15) << meshratio[count] << "\n";
  //   // Close the file
  //   Mesh_file.close();
  // }

  //  else std::cout << "Unable to open file";
}
