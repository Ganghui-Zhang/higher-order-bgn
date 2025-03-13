#include <config.h>

#include <cmath>

#include <dune/common/fvector.hh>
#include <dune/common/parametertree.hh>
#include <dune/common/parametertreeparser.hh>

#include <dune/foamgrid/foamgrid.hh>

#include <dune/gmsh4/gmsh4reader.hh>
#include <dune/gmsh4/gridcreators/discontinuousgridcreator.hh>

#include "eoc_test.hh"
#include <dune/curvedgrid/curvedgrid.hh>
#include <dune/geometry/utility/algorithms.hh>
#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>
#include "dune/foamgrid/foamgrid/foamgridentity.hh"
#include "dune/vtk/datacollectors/lagrangedatacollector.hh"
#include "dune/vtk/types.hh"
#include <dune/curvedgrid/geometries/implicitsurface.hh>
#include <dune/grid/onedgrid.hh>
#include <dune/grid/uggrid.hh>

#include <dune/meshdist/hausdorffdistance.hh>

using namespace Dune;
using namespace std;


template <class T = double>
class EllipsoidProjection
{
  using Domain = FieldVector<T,3>;
  using Jacobian = FieldMatrix<T,3,3>;

  // Implicit function representation of the ellipsoid surface
  struct Phi
  {
    T a_, b_, c_;

    // phi(x,y,z) = (x/a)^2 + (y/b)^2 + (z/c)^2 = 1
    T operator() (const FieldVector<T,3>& x) const
    {
      return x[0]*x[0]/(a_*a_) + x[1]*x[1]/(b_*b_) + x[2]*x[2]/(c_*c_) - 1;
    }

    // grad(phi)
    friend auto derivative (Phi phi)
    {
      return [a=phi.a_,b=phi.b_,c=phi.c_](const Domain& x) -> Domain
      {
        return { T(2*x[0]/(a*a)), T(2*x[1]/(b*b)), T(2*x[2]/(c*c)) };
      };
    }
  };

public:
  /// \brief Constructor of ellipsoid by major axes
  EllipsoidProjection (T a, T b, T c)
    : a_(a)
    , b_(b)
    , c_(c)
    , implicit_(Phi{a,b,c},100)
  {}

  /// \brief project the coordinate to the ellipsoid
  Domain operator() (const Domain& x) const
  {
    return implicit_(x);
  }

private:
  T a_, b_, c_;
  SimpleImplicitSurfaceProjection<Phi> implicit_;
};


int main(int argc, char *argv[])
{

  MPIHelper::instance(argc, argv);

  std::string inifile = "mcf.ini";
  if (argc > 1)
    inifile = argv[1];

  ParameterTree pt;
  ParameterTreeParser::readINITree(inifile, pt);

  //auto hostGrid = Gmsh4Reader<FoamGrid<2,3>>::createGridFromFile( DUNE_GRID_PATH "ellipsoid.msh");
  //auto hostGrid = Gmsh4Reader<FoamGrid<2,3>>::createGridFromFile(DUNE_GRID_PATH "good_ellipsoid_rough.msh");
  //auto hostGrid = Gmsh4Reader<FoamGrid<2,3>>::createGridFromFile(DUNE_GRID_PATH "good_ellipsoid_very_rough.msh");
  auto hostGrid = Gmsh4Reader<FoamGrid<2,3>>::createGridFromFile(DUNE_GRID_PATH "sphere_very_rough.msh");

  hostGrid->globalRefine(3);
  
  double a = 2;
  double b = 1;
  double c = 1;
  EllipsoidProjection Ellipsoid(a, b, c);


 int kg = pt.get<int>("grid.kg", 1);
 int ko = pt.get<int>("grid.ko", 1);

 double tau_ini = pt.get<double>("adapt.tau_ini", 0.05);

 using namespace Dune::Functions::BasisFactory;
 auto feBasis0 = makeBasis(hostGrid->levelGridView(0), 
    composite(
    power<3>(lagrange(kg), flatInterleaved()),
    lagrange(kg),
    flatLexicographic()));
 auto solution0 = run(pt,  tau_ini/std::pow(2, (kg + 1)*0),  feBasis0, Ellipsoid);
 auto positionBasis0 = Dune::Functions::subspaceBasis(feBasis0, Dune::Indices::_0); 
 auto X0 = Dune::Functions::makeDiscreteGlobalBasisFunction<Dune::FieldVector<double,3>>(positionBasis0, solution0); 

auto feBasis1 = makeBasis(hostGrid->levelGridView(1), 
    composite(
    power<3>(lagrange(kg), flatInterleaved()),
    lagrange(kg),
    flatLexicographic()));
auto solution1 = run(pt,  tau_ini/std::pow(2, (kg + 1)*1), feBasis1, Ellipsoid);
auto positionBasis1 = Dune::Functions::subspaceBasis(feBasis1, Dune::Indices::_0); 
auto X1 = Dune::Functions::makeDiscreteGlobalBasisFunction<Dune::FieldVector<double,3>>(positionBasis1, solution1); 

auto error0 = Dune::MeshDist::symmetricMeanSquareError(X1,X0);
std::cout << "error0  = " << error0  << std::endl;

auto feBasis2 = makeBasis(hostGrid->levelGridView(2), 
    composite(
    power<3>(lagrange(kg), flatInterleaved()),
    lagrange(kg),
    flatLexicographic()));
auto solution2 = run(pt,  tau_ini/std::pow(2, (kg + 1)*2), feBasis2, Ellipsoid);
auto positionBasis2 = Dune::Functions::subspaceBasis(feBasis2, Dune::Indices::_0); 
auto X2 = Dune::Functions::makeDiscreteGlobalBasisFunction<Dune::FieldVector<double,3>>(positionBasis2, solution2); 

auto error1 = Dune::MeshDist::symmetricMeanSquareError(X2,X1);
std::cout << "error1  = " << error1  << std::endl;
double order0 = log(error0/error1)/log(2);
std::cout << "order0 = " << order0 << std::endl;


auto feBasis3 = makeBasis(hostGrid->levelGridView(3), 
    composite(
    power<3>(lagrange(kg), flatInterleaved()),
    lagrange(kg),
    flatLexicographic()));
auto solution3 = run(pt,  tau_ini/std::pow(2, (kg + 1)*3), feBasis3, Ellipsoid);
auto positionBasis3 = Dune::Functions::subspaceBasis(feBasis3, Dune::Indices::_0); 
auto X3 = Dune::Functions::makeDiscreteGlobalBasisFunction<Dune::FieldVector<double,3>>(positionBasis3, solution3);

auto error2 = Dune::MeshDist::symmetricMeanSquareError(X3,X2);
std::cout << "error2  = " << error2  << std::endl;
double order1 = log(error1/error2)/log(2);
std::cout << "order1 = " << order1 << std::endl;

}
