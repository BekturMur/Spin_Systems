#pragma once

#include <cstddef>
#include <vector>

namespace spinsim {

/// A Bessel-function sequence J_0(x)..J_{n-1}(x) or I_0(x)..I_{n-1}(x).
struct BesselSequence {
  /// values[k] is the order-k function at the requested argument.
  std::vector<double> values;

  /// Number of leading entries that are non-zero. Beyond this the true values
  /// are below the smallest representable double and are returned as exact
  /// zeros.
  ///
  /// This plays the role of the legacy `NCALC` out-parameter, but with a
  /// meaning that is a property of the answer rather than of the algorithm:
  /// Cody's routines returned NCALC < NB both when values underflowed and when
  /// the routine ran short of precision. Here underflow is simply reported,
  /// and genuine failures throw.
  std::size_t nonzero = 0;
};

/// Bessel functions of the first kind, J_0(x) .. J_{count-1}(x), for x >= 0.
///
/// Replaces rjbesl.f (Cody, ACM TOMS 715), called from chstepsPDDGnmr.f:98 to
/// build the real-time Chebyshev coefficients. Computed by Miller's backward
/// recurrence with the normalisation J_0 + 2(J_2 + J_4 + ...) = 1, which is the
/// same method the Fortran used; agreement with it is asserted by the tests.
///
/// Throws std::invalid_argument for x < 0 or count == 0.
[[nodiscard]] BesselSequence besselJ(double x, std::size_t count);

/// Modified Bessel functions of the first kind, I_0(x) .. I_{count-1}(x), for
/// x >= 0, unscaled (the legacy call sites pass IZE = 1).
///
/// Replaces ribesl.f, called from chimstepsPDDGnmr.f:91 for imaginary-time
/// propagation. Same backward recurrence, normalised by
/// I_0 + 2(I_1 + I_2 + ...) = exp(x).
///
/// Note that unscaled I_n(x) grows like exp(x)/sqrt(2 pi x); for large x the
/// values overflow. Throws std::overflow_error in that case rather than
/// returning infinities.
[[nodiscard]] BesselSequence besselI(double x, std::size_t count);

}  // namespace spinsim
