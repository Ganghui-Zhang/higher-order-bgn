#include <config.h>

#include <cmath>

#include <dune/common/fvector.hh>
#include <dune/common/parametertree.hh>
#include <dune/common/parametertreeparser.hh>

#include <dune/foamgrid/foamgrid.hh>

#include <dune/gmsh4/gmsh4reader.hh>
#include <dune/gmsh4/gridcreators/discontinuousgridcreator.hh>

#include "Curve.hh"
#include <dune/curvedgrid/curvedgrid.hh>
#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>
#include "dune/foamgrid/foamgrid/foamgridentity.hh"
#include "dune/grid/common/rangegenerators.hh"
#include "dune/vtk/datacollectors/lagrangedatacollector.hh"
#include "dune/vtk/types.hh"
#include "implicitsurface.hh"
#include <dune/grid/onedgrid.hh>
#include <dune/grid/uggrid.hh>

using namespace Dune;

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
  /// \brief Constructor of ellipse by major axes
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


class PerturbedCircle
{
public:
  PerturbedCircle(double radius, double twist, double period)
    : radius_(radius)
    , twist_(twist)
    , period_(period)
  {}

  FieldVector<double,2> operator()(FieldVector<double,2> const& x) const
  {
    using std::sqrt;
    using std::atan2;
    double alpha = atan2(x[1], x[0]);
    return x * (radius(alpha)/x.two_norm());
  }

private:
  double radius(double alpha) const
  {
    using std::cos;
    using std::sin;
    double factor = twist_ +  (cos(period_*alpha));
    return radius_ * factor;
  }

private:
  double radius_;
  double twist_;
  double period_;
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
  double mesh_size = 128.0;
  for (std::size_t i = 0; i < mesh_size; ++i){
    factory.insertVertex({a*std::cos(2*M_PI*i/mesh_size), b*std::sin(2*M_PI*i/mesh_size), 0});
  }

  // flower case
  // int k = 6;
  // GridFactory<FoamGrid<1,2>> factory;
  // double mesh_size = 128.0;
  // for (std::size_t i = 0; i < mesh_size; ++i){
  //   factory.insertVertex({(2+std::cos(k*2*M_PI*i/mesh_size))*std::cos(2*M_PI*i/mesh_size), (2+std::cos(k*2*M_PI*i/mesh_size))*std::sin(2*M_PI*i/mesh_size), 0});
  // }

  for (unsigned int i = 0; i < mesh_size; ++i){
    unsigned int next = (i+1) % (unsigned int)(mesh_size);
    factory.insertElement(GeometryTypes::simplex(1),{i, next});
  }

   auto hostGrid = factory.createGrid();

  // Construct a (flat) host grid from a grid file
 // auto hostGrid = Gmsh4Reader<FoamGrid<1,2>>::createGridFromFile(DUNE_GRID_PATH "ellipse.msh");
 
  int refinements = pt.get<int>("grid.refinement_levels", 1);
  hostGrid->globalRefine(refinements);

   //double radius = 1;
   //double twist = 2;
   //int period = 6;
  //  PerturbedCircle flower(radius, twist, period);

  EllipseProjection Ellipse(a, b);

  // double r = 1;
  // Circle circle(r);

  run(pt, *hostGrid, Ellipse);
  //run(pt, *hostGrid, circle);
  //run(pt, *hostGrid, flower);
}
