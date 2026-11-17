// Copyright 2026 Sascha Buettner
#pragma once
#include <array>
#include <cmath>

namespace tam::ocd::point_mass
{
/// Second-order forward-mode number over N independent variables: value, gradient and the
/// Hessian (upper triangle, row-major packed: (0,0) (0,1) .. (0,N-1) (1,1) .. (N-1,N-1)).
///
/// Every elementary operation applies the exact chain rule, so a function templated on the
/// scalar type yields its exact Jacobian and Hessian in one pass; for N = 7 (the five states and
/// two controls of the point mass) one operation costs about 100 flops and the loops are
/// fixed-size, so the compiler unrolls and vectorises them.
template <int N>
struct Dual2
{
  static constexpr int n = N;
  static constexpr int nh = N * (N + 1) / 2;
  double v{0.0};
  std::array<double, N> g{};
  std::array<double, nh> h{};

  Dual2() = default;
  Dual2(double value) : v(value) {}  // NOLINT(runtime/explicit): constants mix freely
  static Dual2 variable(double value, int i)
  {
    Dual2 d(value);
    d.g[i] = 1.0;
    return d;
  }
};

/// Index of (i, j), i <= j, in the packed upper triangle.
constexpr int packed_index(int n, int i, int j) { return i * n - i * (i - 1) / 2 + (j - i); }

/// f(a) with f0 = f(a.v), f1 = f'(a.v), f2 = f''(a.v).
template <int N>
inline Dual2<N> chain(const Dual2<N> & a, double f0, double f1, double f2)
{
  Dual2<N> r;
  r.v = f0;
  for (int i = 0; i < N; ++i) r.g[i] = f1 * a.g[i];
  int k = 0;
  for (int i = 0; i < N; ++i) {
    const double gi = f2 * a.g[i];
    for (int j = i; j < N; ++j, ++k) r.h[k] = f1 * a.h[k] + gi * a.g[j];
  }
  return r;
}

/// f(a, b) with the value f, gradient (fa, fb) and Hessian (faa, fab, fbb) of f at (a.v, b.v).
template <int N>
inline Dual2<N> chain2(
  const Dual2<N> & a, const Dual2<N> & b, double f, double fa, double fb, double faa, double fab,
  double fbb)
{
  Dual2<N> r;
  r.v = f;
  for (int i = 0; i < N; ++i) r.g[i] = fa * a.g[i] + fb * b.g[i];
  int k = 0;
  for (int i = 0; i < N; ++i) {
    const double ai = a.g[i], bi = b.g[i];
    const double c_a = faa * ai + fab * bi, c_b = fab * ai + fbb * bi;
    for (int j = i; j < N; ++j, ++k)
      r.h[k] = fa * a.h[k] + fb * b.h[k] + c_a * a.g[j] + c_b * b.g[j];
  }
  return r;
}

// ---------------------------------------------------------------------------- arithmetic

template <int N>
inline Dual2<N> operator-(const Dual2<N> & a)
{
  Dual2<N> r;
  r.v = -a.v;
  for (int i = 0; i < N; ++i) r.g[i] = -a.g[i];
  for (int k = 0; k < Dual2<N>::nh; ++k) r.h[k] = -a.h[k];
  return r;
}

template <int N>
inline Dual2<N> operator+(const Dual2<N> & a, const Dual2<N> & b)
{
  Dual2<N> r;
  r.v = a.v + b.v;
  for (int i = 0; i < N; ++i) r.g[i] = a.g[i] + b.g[i];
  for (int k = 0; k < Dual2<N>::nh; ++k) r.h[k] = a.h[k] + b.h[k];
  return r;
}

template <int N>
inline Dual2<N> operator-(const Dual2<N> & a, const Dual2<N> & b)
{
  Dual2<N> r;
  r.v = a.v - b.v;
  for (int i = 0; i < N; ++i) r.g[i] = a.g[i] - b.g[i];
  for (int k = 0; k < Dual2<N>::nh; ++k) r.h[k] = a.h[k] - b.h[k];
  return r;
}

template <int N>
inline Dual2<N> operator*(const Dual2<N> & a, const Dual2<N> & b)
{
  Dual2<N> r;
  r.v = a.v * b.v;
  for (int i = 0; i < N; ++i) r.g[i] = a.v * b.g[i] + b.v * a.g[i];
  int k = 0;
  for (int i = 0; i < N; ++i) {
    const double ai = a.g[i], bi = b.g[i];
    for (int j = i; j < N; ++j, ++k)
      r.h[k] = a.v * b.h[k] + b.v * a.h[k] + ai * b.g[j] + bi * a.g[j];
  }
  return r;
}

template <int N>
inline Dual2<N> operator/(const Dual2<N> & a, const Dual2<N> & b)
{
  const double inv = 1.0 / b.v;
  return a * chain(b, inv, -inv * inv, 2.0 * inv * inv * inv);
}

template <int N>
inline Dual2<N> operator+(const Dual2<N> & a, double c)
{
  Dual2<N> r = a;
  r.v += c;
  return r;
}
template <int N>
inline Dual2<N> operator+(double c, const Dual2<N> & a)
{
  return a + c;
}
template <int N>
inline Dual2<N> operator-(const Dual2<N> & a, double c)
{
  Dual2<N> r = a;
  r.v -= c;
  return r;
}
template <int N>
inline Dual2<N> operator-(double c, const Dual2<N> & a)
{
  Dual2<N> r = -a;
  r.v += c;
  return r;
}
template <int N>
inline Dual2<N> operator*(const Dual2<N> & a, double c)
{
  Dual2<N> r;
  r.v = a.v * c;
  for (int i = 0; i < N; ++i) r.g[i] = a.g[i] * c;
  for (int k = 0; k < Dual2<N>::nh; ++k) r.h[k] = a.h[k] * c;
  return r;
}
template <int N>
inline Dual2<N> operator*(double c, const Dual2<N> & a)
{
  return a * c;
}
template <int N>
inline Dual2<N> operator/(const Dual2<N> & a, double c)
{
  return a * (1.0 / c);
}
template <int N>
inline Dual2<N> operator/(double c, const Dual2<N> & a)
{
  const double inv = 1.0 / a.v;
  return chain(a, c * inv, -c * inv * inv, 2.0 * c * inv * inv * inv);
}

template <int N>
inline Dual2<N> & operator+=(Dual2<N> & a, const Dual2<N> & b)
{
  return a = a + b;
}
template <int N>
inline Dual2<N> & operator-=(Dual2<N> & a, const Dual2<N> & b)
{
  return a = a - b;
}
template <int N>
inline Dual2<N> & operator+=(Dual2<N> & a, double c)
{
  a.v += c;
  return a;
}

// ---------------------------------------------------------------------------- functions

template <int N>
inline Dual2<N> sin(const Dual2<N> & a)
{
  const double s = std::sin(a.v), c = std::cos(a.v);
  return chain(a, s, c, -s);
}
template <int N>
inline Dual2<N> cos(const Dual2<N> & a)
{
  const double s = std::sin(a.v), c = std::cos(a.v);
  return chain(a, c, -s, -c);
}
template <int N>
inline Dual2<N> sqrt(const Dual2<N> & a)
{
  const double r = std::sqrt(a.v);
  return chain(a, r, 0.5 / r, -0.25 / (r * a.v));
}
/// a^p for a constant exponent (a > 0 where p is not an integer).
template <int N>
inline Dual2<N> pow(const Dual2<N> & a, double p)
{
  const double f = std::pow(a.v, p);
  const double f1 = p * std::pow(a.v, p - 1.0);
  const double f2 = p * (p - 1.0) * std::pow(a.v, p - 2.0);
  return chain(a, f, f1, f2);
}
template <int N>
inline Dual2<N> fabs(const Dual2<N> & a)
{
  const double s = a.v > 0.0 ? 1.0 : (a.v < 0.0 ? -1.0 : 0.0);
  return chain(a, std::fabs(a.v), s, 0.0);
}
/// max(a, c): the branch of the larger value (CasADi fmax: the derivative of a only where a > c).
template <int N>
inline Dual2<N> fmax(const Dual2<N> & a, double c)
{
  if (a.v > c) return a;
  return Dual2<N>(c);
}
template <int N>
inline Dual2<N> square(const Dual2<N> & a)
{
  return chain(a, a.v * a.v, 2.0 * a.v, 2.0);
}

inline double square(double a) { return a * a; }

/// Value of a scalar or a dual number.
inline double value(double a) { return a; }
template <int N>
inline double value(const Dual2<N> & a)
{
  return a.v;
}

}  // namespace tam::ocd::point_mass
