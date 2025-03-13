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
class EllipseProjection
{
  using Domain = FieldVector<T,2>;
  using Jacobian = FieldMatrix<T,2,2>;


  // Implicit function representation of the ellipsoid surface
  struct Phi
  {
    T a_, b_;

    // phi(x,y) = (x/a)^2 + (y/b)^2 = 1
    T operator() (const FieldVector<T,2>& x) const
    {
      return x[0]*x[0]/(a_*a_) + x[1]*x[1]/(b_*b_) - 1;
    }

    // grad(phi)
    friend auto derivative (Phi phi)
    {
      return [a=phi.a_,b=phi.b_](const Domain& x) -> Domain
      {
        return { T(2*x[0]/(a*a)), T(2*x[1]/(b*b)) };
      };
    }
  };

public:
  /// \brief Constructor of ellipse by major axes
  EllipseProjection (T a, T b)
    : a_(a)
    , b_(b)
    , implicit_(Phi{a,b},100)
  {}

  /// \brief project the coordinate to the ellipsoid
  Domain operator() (const Domain& x) const
  {
    return implicit_(x);
  }

private:
  T a_, b_;
  SimpleImplicitSurfaceProjection<Phi> implicit_;
};



int main(int argc, char *argv[])
{

  MPIHelper::instance(argc, argv);

  std::string inifile = "initial.ini";
  if (argc > 1)
    inifile = argv[1];

  ParameterTree pt;
  ParameterTreeParser::readINITree(inifile, pt);

  // Construct a (flat) host grid by GridFactory
  // ellipse case

  double a = 2;
  double b = 1;

   GridFactory<FoamGrid<1,2>> factory;
  int mesh_size = 32;
  for (std::size_t i = 0; i < mesh_size; ++i){
    factory.insertVertex({a*std::cos(2*M_PI*i/mesh_size), b*std::sin(2*M_PI*i/mesh_size), 0});
  }

  for (unsigned int i = 0; i < mesh_size; ++i){
    unsigned int next = (i+1) % (unsigned int)(mesh_size);
    factory.insertElement(GeometryTypes::simplex(1),{i, next});
  }

auto hostGrid = factory.createGrid();

// Construct a (flat) host grid from a grid file
// auto hostGrid = Gmsh4Reader<FoamGrid<1,2>>::createGridFromFile(DUNE_GRID_PATH "good_ellipse.msh");

hostGrid->globalRefine(6);

EllipseProjection Ellipse(a, b);
 
//std::vector<double> error; 
// std::vector<double> timestep;

// for(int i = 0; i < 5; i++){
//   timestep[i] = 0.005/std::pow(2, 3*i);
// }

 int kg = pt.get<int>("grid.kg", 1);
 int ko = pt.get<int>("grid.ko", 1);

 double tau_ini = pt.get<double>("adapt.tau_ini", 0.05);

 using namespace Dune::Functions::BasisFactory;
 auto feBasis0 = makeBasis(hostGrid->levelGridView(0), 
    composite(
    power<2>(lagrange(kg), flatInterleaved()),
    lagrange(kg),
    flatLexicographic()));
 auto solution0 = run(pt,  tau_ini/std::pow(2, (kg + 1)*0),  feBasis0, Ellipse);
 auto positionBasis0 = Dune::Functions::subspaceBasis(feBasis0, Dune::Indices::_0); 
 auto X0 = Dune::Functions::makeDiscreteGlobalBasisFunction<Dune::FieldVector<double,2>>(positionBasis0, solution0); 

auto feBasis1 = makeBasis(hostGrid->levelGridView(1), 
    composite(
    power<2>(lagrange(kg), flatInterleaved()),
    lagrange(kg),
    flatLexicographic()));
auto solution1 = run(pt,  tau_ini/std::pow(2, (kg + 1)*1), feBasis1, Ellipse);
auto positionBasis1 = Dune::Functions::subspaceBasis(feBasis1, Dune::Indices::_0); 
auto X1 = Dune::Functions::makeDiscreteGlobalBasisFunction<Dune::FieldVector<double,2>>(positionBasis1, solution1); 

//auto error0 = Dune::MeshDist::symmetricMeanSquareError(X1,X0);
auto error0 = Dune::MeshDist::symmetricMeanError(X1,X0);
//auto error0 = Dune::MeshDist::symmetricHausdorffDistance(X1,X0);
std::cout << "error0  = " << error0  << std::endl;

auto feBasis2 = makeBasis(hostGrid->levelGridView(2), 
    composite(
    power<2>(lagrange(kg), flatInterleaved()),
    lagrange(kg),
    flatLexicographic()));
auto solution2 = run(pt,  tau_ini/std::pow(2, (kg + 1)*2), feBasis2, Ellipse);
auto positionBasis2 = Dune::Functions::subspaceBasis(feBasis2, Dune::Indices::_0); 
auto X2 = Dune::Functions::makeDiscreteGlobalBasisFunction<Dune::FieldVector<double,2>>(positionBasis2, solution2); 

//auto error1 = Dune::MeshDist::symmetricMeanSquareError(X2,X1);
auto error1 = Dune::MeshDist::symmetricMeanError(X2,X1);
//auto error1 = Dune::MeshDist::symmetricHausdorffDistance(X2,X1);
std::cout << "error1  = " << error1  << std::endl;
double order0 = log(error0/error1)/log(2);
std::cout << "order0 = " << order0 << std::endl;

auto feBasis3 = makeBasis(hostGrid->levelGridView(3), 
    composite(
    power<2>(lagrange(kg), flatInterleaved()),
    lagrange(kg),
    flatLexicographic()));
auto solution3 = run(pt,  tau_ini/std::pow(2, (kg + 1)*3), feBasis3, Ellipse);
auto positionBasis3 = Dune::Functions::subspaceBasis(feBasis3, Dune::Indices::_0); 
auto X3 = Dune::Functions::makeDiscreteGlobalBasisFunction<Dune::FieldVector<double,2>>(positionBasis3, solution3); 

//auto error2 = Dune::MeshDist::symmetricMeanSquareError(X3,X2);
auto error2 = Dune::MeshDist::symmetricMeanError(X3,X2);
//auto error2 = Dune::MeshDist::symmetricHausdorffDistance(X3,X2);
std::cout << "error2  = " << error2  << std::endl;
double order1 = log(error1/error2)/log(2);
std::cout << "order1 = " << order1 << std::endl;

auto feBasis4 = makeBasis(hostGrid->levelGridView(4), 
    composite(
    power<2>(lagrange(kg), flatInterleaved()),
    lagrange(kg),
    flatLexicographic()));
auto solution4 = run(pt,  tau_ini/std::pow(2, (kg + 1)*4), feBasis4, Ellipse);
auto positionBasis4 = Dune::Functions::subspaceBasis(feBasis4, Dune::Indices::_0); 
auto X4 = Dune::Functions::makeDiscreteGlobalBasisFunction<Dune::FieldVector<double,2>>(positionBasis4, solution4); 

//auto error3 = Dune::MeshDist::symmetricMeanSquareError(X4,X3);
auto error3 = Dune::MeshDist::symmetricMeanError(X4,X3);
//auto error3 = Dune::MeshDist::symmetricHausdorffDistance(X4,X3);
std::cout << "error3  = " << error3  << std::endl;
double order2 = log(error2/error3)/log(2);
std::cout << "order2 = " << order2 << std::endl;

auto feBasis5 = makeBasis(hostGrid->levelGridView(5), 
    composite(
    power<2>(lagrange(kg), flatInterleaved()),
    lagrange(kg),
    flatLexicographic()));
auto solution5 = run(pt,  tau_ini/std::pow(2, (kg + 1)*5), feBasis5, Ellipse);
auto positionBasis5 = Dune::Functions::subspaceBasis(feBasis5, Dune::Indices::_0); 
auto X5 = Dune::Functions::makeDiscreteGlobalBasisFunction<Dune::FieldVector<double,2>>(positionBasis5, solution5); 

//auto error4 = Dune::MeshDist::symmetricMeanSquareError(X5,X4);
auto error4 = Dune::MeshDist::symmetricMeanError(X5,X4);
//auto error4 = Dune::MeshDist::symmetricHausdorffDistance(X5,X4);
std::cout << "error4  = " << error4  << std::endl;
double order3 = log(error3/error4)/log(2);
std::cout << "order3 = " << order3 << std::endl;

}
