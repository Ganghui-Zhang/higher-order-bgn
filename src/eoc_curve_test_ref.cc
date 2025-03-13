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


template <class T = double>
class Circle
{
  using Domain = FieldVector<T,2>;
  using Jacobian = FieldMatrix<T,2,2>;

  // Implicit function representation of the ellipsoid surface
  struct Phi
  {
    T r_;

    // phi(x,y) = (x/r)^2 + (y/r)^2 = 1
    T operator() (const FieldVector<T,2>& x) const
    {
      return x[0]*x[0]/(r_*r_) + x[1]*x[1]/(r_*r_) - 1;
    }

    // grad(phi)
    friend auto derivative (Phi phi)
    {
      return [r=phi.r_](const Domain& x) -> Domain
      {
        return { T(2*x[0]/(r*r)), T(2*x[1]/(r*r)) };
      };
    }
  };

public:
  /// \brief 
  Circle (T r)
    : r_(r)
    , implicit_(Phi{r},100)
  {}

  /// \brief project the coordinate to the ellipsoid
  Domain operator() (const Domain& x) const
  {
    return implicit_(x);
  }

private:
  T r_;
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

  double a = 1;
  double b = 1;

  double r_ref = sqrt(1-2*0.05);

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

int kg = pt.get<int>("grid.kg", 1);
int ko = pt.get<int>("grid.ko", 1);

double tau_ini = pt.get<double>("adapt.tau_ini", 0.05);

double error0 = 0.0;
double error1 = 0.0; double order0; 

using namespace Dune::Functions::BasisFactory;
 auto feBasis0 = makeBasis(hostGrid->levelGridView(0), 
    composite(
    power<2>(lagrange(kg), flatInterleaved()),
    lagrange(kg),
    flatLexicographic()));
 auto solution0 = run(pt,  tau_ini/std::pow(2, (kg + 1)*0),  feBasis0, Ellipse);
 auto positionBasis0 = Dune::Functions::subspaceBasis(feBasis0, Dune::Indices::_0); 
 auto X0 = Dune::Functions::makeDiscreteGlobalBasisFunction<Dune::FieldVector<double,2>>(positionBasis0, solution0);
 auto localX0 = localFunction(X0); 
 
 double max_radius0 = 0;

 // compute error of X_ref and X0
 for (auto const& e0 : elements(feBasis0.gridView()))
 { 
     localX0.bind(e0);
     Dune::LocalFunctionGeometry geometry0{referenceElement(e0), localX0};
    for (auto const& qp : QuadratureRules<double,1>::rule(e0.type(),ko)){
     auto x0 = localX0(qp.position());
     //max_radius0 = std::max(max_radius0,x0.two_norm());
     //std::cout << "max_radius0 = " <<  max_radius0 << std::endl;
     auto dist = (x0-r_ref*x0/x0.two_norm()).two_norm2(); 
     const auto dS = geometry0.integrationElement(qp.position()) * qp.weight();
     error0 += dist * dS;
   }  
 }
std::cout << std::fixed << std::setprecision(15) << "error0 = " << std::sqrt(error0) << std::endl;


auto feBasis1 = makeBasis(hostGrid->levelGridView(1), 
composite(
power<2>(lagrange(kg), flatInterleaved()),
lagrange(kg),
flatLexicographic()));
auto solution1 = run(pt,  tau_ini/std::pow(2, (kg + 1)*1),  feBasis1, Ellipse);
auto positionBasis1 = Dune::Functions::subspaceBasis(feBasis1, Dune::Indices::_0); 
auto X1 = Dune::Functions::makeDiscreteGlobalBasisFunction<Dune::FieldVector<double,2>>(positionBasis1, solution1);
auto localX1 = localFunction(X1); 

double max_radius1 = 0;

for (auto const& e1 : elements(feBasis1.gridView()))
 {    
   localX1.bind(e1);
   Dune::LocalFunctionGeometry geometry1{referenceElement(e1), localX1};
   for (auto const& qp : QuadratureRules<double,1>::rule(e1.type(),ko)){   
     auto x1 = localX1(qp.position());
     //max_radius1 = std::max(max_radius1,x1.two_norm());
     //std::cout << "max_radius1 = " <<  max_radius1 << std::endl;
     auto dist = (x1-r_ref*x1/x1.two_norm()).two_norm2(); 
     const auto dS = geometry1.integrationElement(qp.position()) * qp.weight();
     error1 += dist * dS;
   }  
 }
std::cout << std::fixed << std::setprecision(15) << "error1 = " << std::sqrt(error1) << std::endl;
std::cout << std::fixed << std::setprecision(15) << "Order0 = " << log(error0/error1)/log(2) << std::endl;

}
