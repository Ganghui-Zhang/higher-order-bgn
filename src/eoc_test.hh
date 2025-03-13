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
  LocalAssembler(Dune::ParameterTree const& pt, const double& tau_)
    : //tau_(pt.get<double>("adapt.timestep",0.1))
     quadOrder_(pt.get<int>("solution.quad_order",10))
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
      Dune::FieldVector<double, n.size()> fAtQP = X_e(qp.position());

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


template < class Basis, class InitialSurface>
auto run(Dune::ParameterTree const& pt, const double& tau,  const Basis& feBasis, InitialSurface const& initialSurface)
{
  
  // int kg = pt.get<int>("grid.kg", 1);
  // int ko = pt.get<int>("grid.ko", 10);

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


  Assembler assembler;
  LocalAssembler localAssembler{pt,tau};

  double startTime = pt.get<double>("adapt.start_time", 0.0);
  double endTime = pt.get<double>("adapt.end_time", 0.03);
  //double tau = pt.get<double>("adapt.timestep", 0.1);

  using std::sqrt;
  double tol = sqrt(std::numeric_limits<double>::epsilon());
 

  double t = startTime;

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
  }

  return solution;
}
