#ifndef DUNE_CURVEDGRID_IMPLICIT_SURFACE_PROJECTION_HH
#define DUNE_CURVEDGRID_IMPLICIT_SURFACE_PROJECTION_HH

#include <cmath>
#include <type_traits>

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/common/parametertree.hh>
#include <dune/common/typeutilities.hh>
#include <dune/curvedgrid/gridfunctions/analyticgridfunction.hh>

namespace Dune {

/// \brief Closest-point projection to surface given by implicit function
/**
 * Surface S is given by zero-levelset of function F.
 * We assume that F is differentiable in order to evaluate normals and to
 * do an iterative projection.
 **/
template <class F>
class ImplicitSurfaceProjection
{
  using Functor = F;
  using DFunctor = std::decay_t<decltype(derivative(std::declval<F>()))>;

public:
  /// \brief Constructor for a given differentiable function with surface = { x : phi(x) == 0 }
  ImplicitSurfaceProjection (const F& phi, int maxIter = 10)
    : maxIter_(maxIter)
    , phi_(phi)
    , dphi_(derivative(phi_))
  {}

  /// \brief Evaluation of the closest-point projection
  /**
   * Iterative algorithm proposed by Demlow and Dziuk in
   * AN ADAPTIVE FINITE ELEMENT METHOD FOR THE LAPLACE-BELTRAMI OPERATOR ON IMPLICITLY DEFINED SURFACES
   **/
  template <class Domain>
  Domain operator() (const Domain& x0) const
  {
    using std::sqrt;
    using T = typename Domain::value_type;

    Domain x = x0;
    T tol = 10*std::numeric_limits<T>::epsilon();

    auto phi = phi_(x);
    auto sign = phi > 0 ? 1 : phi < 0 ? -1 : 0;
    auto grad_phi = dphi_(x);
    auto grad_phi_norm2 = grad_phi.two_norm2();

    T err1 = 1, err2 = 1;
    for (int i = 0; i < maxIter_; ++i) {
      auto y = x - grad_phi * (phi/grad_phi_norm2);

      auto dist = sign * (y - x0).two_norm();
      auto grad_phi_y = dphi_(y);

      auto& dx = grad_phi_y;
      dx *= -dist/grad_phi_y.two_norm();
      x = x0 + dx;

      phi = phi_(x);
      grad_phi = dphi_(x);
      grad_phi_norm2 = grad_phi.two_norm2();

      err1 = sqrt(phi*phi / grad_phi_norm2);
      err2 = std::min(
        (grad_phi / sqrt(grad_phi_norm2) - dx / dx.two_norm()).two_norm(),
        (grad_phi / sqrt(grad_phi_norm2) + dx / dx.two_norm()).two_norm() );
      if (err1 + err2 < tol)
        break;
    }

    return x;
  }

  /// \brief Normal vector given by grad(F)/|grad(F)|
  template <class Domain>
  Domain normal (const Domain& x0) const
  {
    auto grad_phi = dphi_(x0);
    return grad_phi / grad_phi.two_norm();
  }

private:
  int maxIter_;
  Functor phi_;
  DFunctor dphi_;
};


/// \brief Closest-point projection to surface given by implicit function using a simple projection algorithm
/**
 * Surface S is given by zero-levelset of function F.
 * We assume that F is differentiable in order to evaluate normals and to
 * do an iterative projection.
 **/
template <class F>
class SimpleImplicitSurfaceProjection
{
  using Functor = F;
  using DFunctor = std::decay_t<decltype(derivative(std::declval<F>()))>;

public:
  /// \brief Constructor for a given differentiable function with surface = { x : phi(x) == 0 }
  SimpleImplicitSurfaceProjection (const F& phi, int maxIter = 10)
    : maxIter_(maxIter)
    , phi_(phi)
    , dphi_(derivative(phi_))
  {}

  /// \brief Evaluation of the closest-point projection
  /**
   * Simple iterative algorithm proposed by Ingo Nitschke in
   * Diskretes Äußeres Kalkül (DEC) auf Oberflächen ohne Rand
   **/
  template <class Domain>
  Domain operator() (Domain x) const
  {
    using std::sqrt;
    using T = typename Domain::value_type;

    T tol = 10*std::numeric_limits<T>::epsilon();

    for (int i = 0; i < maxIter_; ++i) {
      auto phi_x = phi_(x);
      auto grad_phi_x = dphi_(x);
      auto normalize = phi_x / grad_phi_x.two_norm2();

      if (sqrt(phi_x * normalize) < tol)
        break;

      x -= grad_phi_x * normalize;
    }

    return x;
  }

  /// \brief Normal vector given by grad(F)/|grad(F)|
  template <class Domain>
  Domain normal (const Domain& x0) const
  {
    auto grad_phi = dphi_(x0);
    return grad_phi / grad_phi.two_norm();
  }

private:
  int maxIter_;
  Functor phi_;
  DFunctor dphi_;
};



/// \brief Closest-point projection to surface given by implicit function
/**
 * Surface S is given by zero-levelset of function F.
 * We assume that F is differentiable in order to evaluate normals and to
 * do an iterative projection.
 **/
template <class Phi, class T = double>
class NewtonImplicitSurfaceProjection
{
  using Functor = Phi;
  using DFunctor = std::decay_t<decltype(derivative(std::declval<Phi>()))>;
  using D2Functor = std::decay_t<decltype(derivative(derivative(std::declval<Phi>())))>;

public:
  /// \brief Constructor for a given differentiable function with surface = { x : phi(x) == 0 }
  NewtonImplicitSurfaceProjection (const Phi& phi, int maxIter = 10)
    : maxIter_(maxIter)
    , phi_(phi)
    , dphi_(derivative(phi_))
    , d2phi_(derivative(dphi_))
  {}

  /// \brief Evaluation of the closest-point projection
  /**
   * Iterative algorithm proposed by Demlow and Dziuk in
   * AN ADAPTIVE FINITE ELEMENT METHOD FOR THE LAPLACE-BELTRAMI OPERATOR ON IMPLICITLY DEFINED SURFACES
   **/
  template <class Domain>
  Domain operator() (const Domain& x0_) const
  {
    using std::sqrt; using std::min;

    T tol = 10*std::numeric_limits<typename Domain::value_type>::epsilon();

    FieldVector<T,3> x0 = x0_;

    auto F = [&x0](auto const& phi_x, auto const& grad_phi_x, auto const& x, auto const& lambda) -> FieldVector<T,4>
    {
      return {
        2*(x[0] - x0[0]) + lambda*grad_phi_x[0],
        2*(x[1] - x0[1]) + lambda*grad_phi_x[1],
        2*(x[2] - x0[2]) + lambda*grad_phi_x[2],
        phi_x
      };
    };

    auto J = [](auto const& grad_phi_x, auto const& H_phi_x, auto const& lambda) -> FieldMatrix<T,4,4>
    {
      return {
        {lambda*H_phi_x[0][0]+2, lambda*H_phi_x[0][1],   lambda*H_phi_x[0][2],   grad_phi_x[0]},
        {lambda*H_phi_x[1][0],   lambda*H_phi_x[1][1]+2, lambda*H_phi_x[1][2],   grad_phi_x[1]},
        {lambda*H_phi_x[2][0],   lambda*H_phi_x[2][1],   lambda*H_phi_x[2][2]+2, grad_phi_x[2]},
        {grad_phi_x[0],          grad_phi_x[1],          grad_phi_x[2],          0}
      };
    };

    FieldVector<T,3> x = x0;
    x -= dphi_(x) * phi_(x) / dphi_(x).two_norm2();

    auto phi = phi_(x);
    auto grad_phi = dphi_(x);
    auto grad_phi_norm2 = grad_phi.two_norm2();

    T lambda = 2*phi/grad_phi_norm2;

    T err1 = 1, err2 = 1;
    for (int i = 0; i < maxIter_; ++i) {
      auto dx = x - x0;

      err1 = sqrt(phi*phi / grad_phi_norm2);
      err2 = min(
        (grad_phi / sqrt(grad_phi_norm2) - dx / dx.two_norm()).two_norm(),
        (grad_phi / sqrt(grad_phi_norm2) + dx / dx.two_norm()).two_norm() );
      if (err1 + err2 < tol)
        break;

      FieldVector<T,4> d;
      J(grad_phi, d2phi_(x), lambda).solve(d, -F(phi, grad_phi, x, lambda));

      x[0] += d[0];
      x[1] += d[1];
      x[2] += d[2];

      lambda += d[3];

      phi = phi_(x);
      grad_phi = dphi_(x);
      grad_phi_norm2 = grad_phi.two_norm2();
    }

    return x;
  }

  /// \brief Normal vector given by grad(F)/|grad(F)|
  template <class Domain>
  Domain normal (const Domain& x0) const
  {
    auto grad_phi = dphi_(x0);
    return grad_phi / grad_phi.two_norm();
  }

private:
  int maxIter_;
  Functor phi_;
  DFunctor dphi_;
  D2Functor d2phi_;
};


/// \brief Construct a grid function representing a projection to an implicit surface
template <class Grid, class F>
auto implicitSurfaceGridFunction (const F& phi, int maxIter = 10)
{
  return analyticGridFunction<Grid>(ImplicitSurfaceProjection<F>{phi, maxIter});
}

/// \brief Construct a grid function representing a projection to an implicit surface by
/// additionally giving the projection implementation as template
template <class Grid, template <class> class Projection, class F>
auto implicitSurfaceGridFunction (const F& phi, int maxIter = 10)
{
  return analyticGridFunction<Grid>(Projection<F>{phi, maxIter});
}

} // end namespace Dune

#endif // DUNE_CURVEDGRID_IMPLICIT_SURFACE_PROJECTION_HH
