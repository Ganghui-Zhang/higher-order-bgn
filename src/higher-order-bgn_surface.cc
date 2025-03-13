#include <config.h>

#include <cmath>

#include <dune/common/fvector.hh>
#include <dune/common/parametertree.hh>
#include <dune/common/parametertreeparser.hh>

#include <dune/foamgrid/foamgrid.hh>

#include <dune/gmsh4/gmsh4reader.hh>

#include "Surface.hh"
#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>
#include "implicitsurface.hh"

using namespace Dune;

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

  /// \brief Normal vector
  Domain normal (const Domain& X) const
  {
    using std::sqrt;
    T x = X[0], y = X[1], z = X[2];
    T a2 = a_*a_, b2 = b_*b_, c2 = c_*c_;

    auto div = sqrt(b2*b2*c2*c2*x*x + a2*a2*c2*c2*y*y + a2*a2*b2*b2*z*z);
    return { b2*c2*x/div , a2*c2*y/div , a2*b2*z/div };
  }

  /// \brief Mean curvature
  T mean_curvature (const Domain& X) const
  {
    using std::sqrt; using std::abs;
    T x = X[0], y = X[1], z = X[2];
    T x2 = x*x, y2 = y*y, z2 = z*z;
    T a2 = a_*a_, b2 = b_*b_, c2 = c_*c_;

    auto div = 2*a2*b2*c2 * power(sqrt(x2/(a2*a2) + y2/(b2*b2) + z2/(c2*c2)), 3);
    return abs(x2 + y2 + z2 - a2 - b2 - c2)/div;
  }

  /// \brief Gaussian curvature
  T gauss_curvature (const Domain& X) const
  {
    T x = X[0], y = X[1], z = X[2];
    T x2 = x*x, y2 = y*y, z2 = z*z;
    T a2 = a_*a_, b2 = b_*b_, c2 = c_*c_;

    auto div = a2*b2*c2 * power(x2/(a2*a2) + y2/(b2*b2) + z2/(c2*c2), 2);
    return T(1)/div;
  }

private:
  T a_, b_, c_;
  SimpleImplicitSurfaceProjection<Phi> implicit_;
};

template <class T = double>
class SphereProjection
{
  using Domain = FieldVector<T,3>;
  using Jacobian = FieldMatrix<T,3,3>;

  // Implicit function representation of the ellipsoid surface
  struct Phi
  {
    T R_;

    // phi(x,y,z) = (x/R)^2 + (y/R)^2 + (z/R)^2 = 1
    T operator() (const FieldVector<T,3>& x) const
    {
      return x[0]*x[0]/(R_*R_) + x[1]*x[1]/(R_*R_) + x[2]*x[2]/(R_*R_) - 1;
    }

    // grad(phi)
    friend auto derivative (Phi phi)
    {
      return [R=phi.R_](const Domain& x) -> Domain
      {
        return { T(2*x[0]/(R*R)), T(2*x[1]/(R*R)), T(2*x[2]/(R*R)) };
      };
    }
  };

public:
  /// \brief Constructor of ellipsoid by major axes
  SphereProjection (T R)
    : R_(R)
    , implicit_(Phi{R},100)
  {}

  /// \brief project the coordinate to the ellipsoid
  Domain operator() (const Domain& x) const
  {
    return implicit_(x);
  }

private:
  T R_;
  SimpleImplicitSurfaceProjection<Phi> implicit_;
};


template <class T = double>
class TorusProjection
{
  using Domain = FieldVector<T,3>;
  using Jacobian = FieldMatrix<T,3,3>;

  // Implicit function representation of the torus surface
  struct Phi
  {
    T R, r;

    // phi(x,y,z) = (sqrt(x^2 + y^2) - R)^2 + z^2 = r^2
    T operator() (const Domain& x) const
    {
      using std::sqrt;
      T phi0 = sqrt(x[0]*x[0] + x[1]*x[1]);
      return (phi0 - R)*(phi0 - R) + x[2]*x[2] - r*r;
    }

    // grad(phi)
    friend auto derivative (const Phi& phi)
    {
      return [R=phi.R,r=phi.r](const Domain& x) -> Domain
      {
        using std::sqrt;
        T phi0 = sqrt(x[0]*x[0] + x[1]*x[1]);
        return {
          -2*R*x[0]/phi0 + 2*x[0] ,
          -2*R*x[1]/phi0 + 2*x[1] ,
            2*x[2]
        };
      };
    }
  };

public:
  // \brief Construction of the torus function by major and minor radius
  TorusProjection (T R, T r)
    : R_(R)
    , r_(r)
  {}

  /// \brief Closest point projection
  Domain operator() (const Domain& x) const
  {
    using std::sqrt;
    auto scale1 = R_ / sqrt(x[0]*x[0] + x[1]*x[1]);
    Domain center{x[0] * scale1, x[1] * scale1, 0};
    Domain out{x[0] - center[0], x[1] - center[1], x[2]};
    out *= r_ / out.two_norm();

    return out + center;
  }


  /// \brief Normal vector
  Domain normal (const Domain& x) const
  {
    using std::sqrt;
    auto X = (*this)(x);
    T factor = (1 - R_/sqrt(X[0]*X[0] + X[1]*X[1]))/r_;
    return { X[0]*factor, X[1]*factor, X[2]/r_ };
  }

  /// \brief Mean curvature
  T mean_curvature (const FieldVector<T,3>& x) const
  {
    using std::sqrt;
    auto X = (*this)(x);
    return (2 - R_/sqrt(X[0]*X[0] + X[1]*X[1]))/(2*r_);
  }

  /// \brief surface area of the torus = 4*pi^2*R*r
  T area () const
  {
    return 4*M_PI*M_PI*R_*r_;
  }

private:
  T R_, r_;
};



int main(int argc, char *argv[])
{
  MPIHelper::instance(argc, argv);

  std::string inifile = "initial.ini";
  if (argc > 1)
    inifile = argv[1];

  ParameterTree pt;
  ParameterTreeParser::readINITree(inifile, pt);

  // Construct a (flat) host grid from a grid file
  //auto hostGrid = Gmsh4Reader<FoamGrid<2,3>>::createGridFromFile(DUNE_GRID_PATH "sphere.msh");
  auto hostGrid = Gmsh4Reader<FoamGrid<2,3>>::createGridFromFile(DUNE_GRID_PATH "ellipsoid.msh");
  //auto hostGrid = Gmsh4Reader<FoamGrid<2,3>>::createGridFromFile(DUNE_GRID_PATH "torus.msh");
  //auto hostGrid = Gmsh4Reader<FoamGrid<2,3>>::createGridFromFile(DUNE_GRID_PATH "sphere_very_rough.msh");

  int refinements = pt.get<int>("grid.refinement_levels", 1);
  hostGrid->globalRefine(refinements);

  // initial surface parametrization
   double a = pt.get<double>("grid.initial.a", 2);
   double b = pt.get<double>("grid.initial.b", 1);
   double c = pt.get<double>("grid.initial.c", 1);
  EllipsoidProjection Ellipsoid(a, b, c);

  //double R = 1;
  //SphereProjection Sphere(R);
  // TorusProjection Torus(R,r);

  run(pt, *hostGrid, Ellipsoid);
  //run(pt, *hostGrid, Sphere);
  //run(pt, *hostGrid, Torus);
}
