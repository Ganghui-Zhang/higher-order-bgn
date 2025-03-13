#pragma once
#include <cmath>
#include <type_traits>
#include <vector>
#include <iostream>
#include <fstream>
using namespace std;

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/common/parametertree.hh>

#include <dune/curvedgrid/curvedgrid.hh>
#include <dune/curvedgrid/gridfunctions/discretegridviewfunction.hh>

#include <dune/curvedgeometry/parametrizedgeometry.hh>

#include <dune/functions/gridfunctions/composedgridfunction.hh>

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

class LocalAssemblerPicard
{
public:
  LocalAssemblerPicard(Dune::ParameterTree const& pt)
    : tau_(pt.get<double>("adapt_sp.timestep",0.1))
    , quadOrder_(pt.get<int>("solution.quad_order",4))
  {}

  static auto perp (Dune::FieldMatrix<double,2,3> const& J){
   return Dune::FieldVector<double,3>{J[0][1] * J[1][2] - J[0][2] * J[1][1],
      J[0][2] * J[1][0] - J[0][0] * J[1][2],
      J[0][0] * J[1][1] - J[0][1] * J[1][0]};
  }
  
  // Compute the stiffness matrix for a single element
  template <class LocalView, class LocalParametrization, class LocalParametrizationIter, class LocalParametrizationMid,
            class Geometry, class GeometryIter, class GeometryMid, class MatrixType, class VectorType>
  void operator()(LocalView const& localView, LocalParametrization const& X_e, LocalParametrizationIter const& Xiter_e, LocalParametrizationMid const& Xmid_e,
                  Geometry const& geometry, GeometryIter const& geometry_iter, GeometryMid const& geometry_mid, MatrixType& elementMatrix, VectorType& elementVector) const
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
    const auto& quad = Dune::QuadratureRules<double,2>::rule(element.type(), quadOrder_);

    // Loop over all quadrature points
    for (const auto& qp : quad)
    {
      // integration element dG_{s-1}
      const auto dS = geometry.integrationElement(qp.position()) * qp.weight();

      // the inverse of the transposed geometry Jacobian
      const auto& Jtinv = geometry.jacobianInverseTransposed(qp.position());
  

      auto const J1 = geometry.jacobianTransposed(qp.position());
      auto const J2 = geometry_mid.jacobianTransposed(qp.position());
      auto const J3 = geometry_iter.jacobianTransposed(qp.position());

      auto n_Picard = (perp(J1) + 4* perp(J2) + perp(J3)) * qp.weight() / (6* dS);

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

          // tau * (grad(v), grad(w)) for sd
           elementMatrix[row][col] += -tau_ * laplace;

          for (int i = 0; i < n_Picard.size(); ++i) {
            const auto mass_normal = shapeValues[j1] * shapeValues[j2] * n_Picard[i] * dS;
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
      Dune::FieldVector<double,3> fAtQP = X_e(qp.position());

      for (std::size_t j1 = 0; j1 < s0; ++j1) {
        for (int i = 0; i < 3; ++i) {
          std::size_t row = node.child(_1).localIndex(j1);

          // (x, normal*w)
          elementVector[row] += fAtQP[i] * shapeValues[j1] * n_Picard[i] * dS;
        } // end i
      } // end j1
    }
  }

private:
  double tau_;
  int quadOrder_;
};


struct AssemblerPicard
{
  // assembler matrix and vector on the grid
  template <class Basis, class Parametrization,  class ParametrizationIter,  class ParametrizationMid, class LocalAssemblerPicard, class BM, class BV>
  void operator()(Basis const& basis, Parametrization const& X, ParametrizationIter const& Xiter, ParametrizationMid const& X_mid, LocalAssemblerPicard const& localAssemblerPicard,
                  Dune::BCRSMatrix<BM>& matrix, Dune::BlockVector<BV>& rhs) const
  {
    // resize matrix and vector and initialize nonzero structure
    this->initialize(basis, matrix);
    this->initialize(basis, rhs);

    Dune::Matrix<double> elementMatrix;
    Dune::BlockVector<double> elementVector;

    auto X_e = localFunction(X);
    auto Xiter_e = localFunction(Xiter);
    auto X_mid_e = localFunction(X_mid);
    auto localView = basis.localView();
  

    for(const auto& e : elements(basis.gridView()))
    {
      X_e.bind(e);
      Xiter_e.bind(e);
      X_mid_e.bind(e);
      localView.bind(e);

      auto const& localFEmid = localView.tree().child(Dune::Indices::_0).child(0).finiteElement();

      // construct the curved geometry from the parametrization
      Dune::LocalFunctionGeometry curvedGeometry{referenceElement(e), X_e};
      Dune::LocalFunctionGeometry curvedGeometryIter{referenceElement(e), Xiter_e};
      Dune::ParametrizedGeometry curvedGeometryMid{referenceElement(e), localFEmid, X_mid_e};

    
      localAssemblerPicard(localView, X_e, Xiter_e, X_mid_e, curvedGeometry, curvedGeometryIter, curvedGeometryMid, elementMatrix, elementVector);
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

/// Compute geometric quantities: surface area and volumn
  template <class Basis, class Parametrization>
  auto SurfaceAreaAndVolume(Dune::ParameterTree const& pt, Basis const& basis, Parametrization const& X)
  {
    double SurfaceArea = 0;
    double Volume = 0;

    auto X_e = localFunction(X);
    auto localView = basis.localView();

    int ko = pt.get<int>("grid.ko", 1);

    for (auto const& e : elements(basis.gridView()))
    {
      X_e.bind(e);
      localView.bind(e);

      // construct the curved geometry from the parametrization
      Dune::LocalFunctionGeometry geometry{referenceElement(e), X_e};

      auto const& quadRule = Dune::QuadratureRules<double,2>::rule(e.type(), ko);
      for (auto const& qp : quadRule)
       {
         const auto dS = geometry.integrationElement(qp.position()) * qp.weight();
         auto const& n = geometry.normal(qp.position());

        //  surface_area = int_{gridview} dx
         SurfaceArea += dS;
          
         Dune::FieldVector<double,3> fAtQP = X_e(qp.position());
         for (int i = 0; i < 3; ++i) {
          //  volume = int_{gridview} 1/3*X_e*normal dx
           Volume += 1.0/3 * fAtQP[i] * n[i] * dS;
         }
      }
    }

    return std::tuple{SurfaceArea,Volume};
  }



template <class HostGrid, class InitialSurface>
void run(Dune::ParameterTree const& pt, HostGrid& hostGrid, InitialSurface const& initialSurface)
{
  int kg = pt.get<int>("grid.kg", 1);

  using namespace Dune::Functions::BasisFactory;
  auto feBasis = makeBasis(hostGrid.leafGridView(), 
    composite(
      power<3>(lagrange(kg), flatInterleaved()),
      lagrange(kg),
      flatLexicographic()));

  auto positionBasis = Dune::Functions::subspaceBasis(feBasis, Dune::Indices::_0);
  auto curvatureBasis = Dune::Functions::subspaceBasis(feBasis, Dune::Indices::_1);
 
  using Matrix = Dune::BCRSMatrix<double>;
  Matrix matrix;
 
  using Vector = Dune::BlockVector<double>;
  Vector rhs;

  // interpolate initial surface to solution vector
  Vector solution, solution_iter;

  auto Xh = Dune::Functions::makeDiscreteGlobalBasisFunction<Dune::FieldVector<double,3>>(positionBasis, solution);
  auto Xh_iter = Dune::Functions::makeDiscreteGlobalBasisFunction<Dune::FieldVector<double,3>>(positionBasis, solution_iter);

  auto Xh_mid = Dune::Functions::makeComposedGridFunction([](auto X_x,auto X_iter_x){
    return (X_x + X_iter_x)/2.0;}, Xh, Xh_iter);

  auto Hh = Dune::Functions::makeDiscreteGlobalBasisFunction<double>(curvatureBasis, solution);

  Dune::Functions::interpolate(positionBasis, solution, initialSurface);

  // the computational domain
  Dune::CurvedGrid grid{hostGrid, Xh};
  auto gridView = grid.leafGridView();

  // create a grid filewriter function to write solution of transformed grid
  // auto dataCollector = Dune::Vtk::LagrangeDataCollector{gridView, kg};
  // using Writer = decltype(Dune::Vtk::UnstructuredGridWriter{dataCollector});
  // auto pvdWriter = Dune::Vtk::PvdWriter<Writer>{dataCollector};
  // auto writeFiles = [&,filename=pt.get<std::string>("grid.output_filename", "grid")](double time)
  // {
  //   pvdWriter.writeTimestep(time, filename + ".vtu");
  // };

  AssemblerPicard assemblerPicard;
  LocalAssemblerPicard localAssemblerPicard{pt};

  double startTime = pt.get<double>("adapt_sp.start_time", 0.0);
  double endTime = pt.get<double>("adapt_sp.end_time", 0.03);
  double tau = pt.get<double>("adapt_sp.timestep", 0.1);

  std::vector<double> time;
  std::vector<double> surface_area;
  std::vector<double> volume;
  std::vector<double> relative_volume_loss;
  std::vector<double> iteration_error_vector;
  std::vector<double> iteration_number_vector;
  std::vector<double> normalized_surface_area;

  using std::sqrt;
  double tol = sqrt(std::numeric_limits<double>::epsilon());

  double t = startTime;
  int total_step = (endTime - startTime)/tau;
  //int maxIter = 100;
  //double iter_tol = 1e-12;
  double iter_tol = (pt.get<double>("adapt_sp.iter_tol",1e-12));
  int maxIter = (pt.get<double>("adapt_sp.maxIter",100));
  int Picard_iteration_number = 0;
  //writeFiles(t);

  // auto[S0,V0] = SurfaceAreaAndVolume(pt,feBasis,Xh);

  // time.push_back(t);
  // surface_area.push_back(S0);
  // normalized_surface_area.push_back(S0/S0);
  // volume.push_back(V0);
  // relative_volume_loss.push_back((V0 - V0)/V0);

  solution_iter = solution;
  while (t + tau < endTime + tol)
  { 
    t += tau;
    //std::cout << std::endl;
    std::cout << "time = " << std::fixed << std::setprecision(15) << t << ", timestep = " << tau << std::endl;
    
     double iter_error = 1e10;

     for (int iter_number = 0; iter_number < maxIter && iter_error > iter_tol; ++iter_number){
        auto solution_iter_last = solution_iter;

        // assemble matrix and vector
        assemblerPicard(feBasis, Xh, Xh_iter, Xh_mid, localAssemblerPicard, matrix, rhs);

      // create linear solver
      auto solver = Dune::UMFPack<Matrix>(matrix, pt.get<int>("solver.verbose",0));

       // solve the linear system
       Dune::InverseOperatorResult statistics;
       solver.apply(solution_iter, rhs, statistics);
   
        // compute iteration error
        solution_iter_last -= solution_iter;
        iter_error = solution_iter_last.two_norm();
        iteration_error_vector.push_back(iter_error);

        //std::cout << std::scientific << std::setprecision(15) << "iteration_error = " << iter_error << std::endl;
        //std::cout << std::fixed << std::setprecision(15) << "iteration_number = " << iter_number << std::endl;

        Picard_iteration_number = iter_number;
     }  
    solution = solution_iter;

    //  auto[S,V] = SurfaceAreaAndVolume(pt,feBasis,Xh);

    //  time.push_back(t);
    //  surface_area.push_back(S);
    //  normalized_surface_area.push_back(S/S0);
    //  volume.push_back(V);
    //  relative_volume_loss.push_back((V - V0)/V0);
    //  iteration_number_vector.push_back(Picard_iteration_number);

    // std::cout << std::endl;
    // std::cout << "iteration = " << Picard_iteration_number << std::endl;

    // write solution to file
    //writeFiles(t);
  }


//  std::ofstream time_file ("Time.txt");
//   if (time_file.is_open())
//   {
//     for (int count = 0; count < total_step; count++)
//         time_file  
//         << std::fixed << std::setprecision(15) << time[count] << "\n";
//     // Close the file
//     time_file.close();
//   }


  // std::ofstream Surface_Area_file ("Normalized_Surface_Area.txt");
  // if (Surface_Area_file.is_open())
  // {
  //   for (int count = 0; count < total_step; count++)
  //       Surface_Area_file  
  //       << std::fixed << std::setprecision(15) << normalized_surface_area[count] << "\n";
  //   // Close the file
  //   Surface_Area_file.close();
  // }
    
  // std::ofstream Volume_file ("Relative_volume.txt");
  // if (Volume_file.is_open())
  // {
  //   for (int count = 0; count < total_step; count++)
  //       Volume_file  
  //       << std::fixed << std::setprecision(15) << relative_volume_loss[count] << "\n";
  //   // Close the file
  //   Volume_file.close();
  // }

  // std::ofstream Iteration_file ("Iteration_number.txt");
  // if (Iteration_file.is_open())
  // {
  //   for (int count = 0; count < total_step; count++)
  //       Iteration_file  
  //       << std::fixed << std::setprecision(15) << iteration_number_vector[count] << "\n";
  //   // Close the file
  //   Iteration_file.close();
  // }

}
