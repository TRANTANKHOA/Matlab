#ifndef STAN_MATH_PRIM_SCAL_PROB_CAUCHY_CDF_HPP
#define STAN_MATH_PRIM_SCAL_PROB_CAUCHY_CDF_HPP

#include <boost/random/cauchy_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <stan/math/prim/scal/meta/OperandsAndPartials.hpp>
#include <stan/math/prim/scal/err/check_consistent_sizes.hpp>
#include <stan/math/prim/scal/err/check_finite.hpp>
#include <stan/math/prim/scal/err/check_not_nan.hpp>
#include <stan/math/prim/scal/err/check_positive_finite.hpp>
#include <stan/math/prim/scal/fun/constants.hpp>
#include <stan/math/prim/scal/fun/log1p.hpp>
#include <stan/math/prim/scal/fun/square.hpp>
#include <stan/math/prim/scal/fun/value_of.hpp>
#include <stan/math/prim/scal/meta/include_summand.hpp>
#include <limits>

namespace stan {

  namespace math {

    /**
     * Calculates the cauchy cumulative distribution function for
     * the given variate, location, and scale.
     *
     * \f$\frac{1}{\pi}\arctan\left(\frac{y-\mu}{\sigma}\right) + \frac{1}{2}\f$
     *
     * @param y A scalar variate.
     * @param mu The location parameter.
     * @param sigma The scale parameter.
     *
     * @return
     */
    template <typename T_y, typename T_loc, typename T_scale>
    typename return_type<T_y, T_loc, T_scale>::type
    cauchy_cdf(const T_y& y, const T_loc& mu, const T_scale& sigma) {
      typedef typename stan::partials_return_type<T_y, T_loc, T_scale>::type
        T_partials_return;

      // Size checks
      if ( !( stan::length(y) && stan::length(mu)
              && stan::length(sigma) ) )
        return 1.0;

      static const char* function("stan::math::cauchy_cdf");

      using stan::math::check_positive_finite;
      using stan::math::check_finite;
      using stan::math::check_not_nan;
      using stan::math::check_consistent_sizes;
      using boost::math::tools::promote_args;
      using stan::math::value_of;

      T_partials_return P(1.0);

      check_not_nan(function, "Random variable", y);
      check_finite(function, "Location parameter", mu);
      check_positive_finite(function, "Scale parameter", sigma);
      check_consistent_sizes(function,
                             "Random variable", y,
                             "Location parameter", mu,
                             "Scale Parameter", sigma);

      // Wrap arguments in vectors
      VectorView<const T_y> y_vec(y);
      VectorView<const T_loc> mu_vec(mu);
      VectorView<const T_scale> sigma_vec(sigma);
      size_t N = max_size(y, mu, sigma);

      OperandsAndPartials<T_y, T_loc, T_scale>
        operands_and_partials(y, mu, sigma);

      // Explicit return for extreme values
      // The gradients are technically ill-defined, but treated as zero
      for (size_t i = 0; i < stan::length(y); i++) {
        if (value_of(y_vec[i]) == -std::numeric_limits<double>::infinity())
          return operands_and_partials.to_var(0.0, y, mu, sigma);
      }

      // Compute CDF and its gradients
      using std::atan;
      using stan::math::pi;

      // Compute vectorized CDF and gradient
      for (size_t n = 0; n < N; n++) {
        // Explicit results for extreme values
        // The gradients are technically ill-defined, but treated as zero
        if (value_of(y_vec[n]) == std::numeric_limits<double>::infinity()) {
          continue;
        }

        // Pull out values
        const T_partials_return y_dbl = value_of(y_vec[n]);
        const T_partials_return mu_dbl = value_of(mu_vec[n]);
        const T_partials_return sigma_inv_dbl = 1.0 / value_of(sigma_vec[n]);

        const T_partials_return z = (y_dbl - mu_dbl) * sigma_inv_dbl;

        // Compute
        const T_partials_return Pn = atan(z) / pi() + 0.5;

        P *= Pn;

        if (!is_constant_struct<T_y>::value)
          operands_and_partials.d_x1[n]
            += sigma_inv_dbl / (pi() * (1.0 + z * z) * Pn);
        if (!is_constant_struct<T_loc>::value)
          operands_and_partials.d_x2[n]
            += - sigma_inv_dbl / (pi() * (1.0 + z * z) * Pn);
        if (!is_constant_struct<T_scale>::value)
          operands_and_partials.d_x3[n]
            += - z * sigma_inv_dbl / (pi() * (1.0 + z * z) * Pn);
      }

      if (!is_constant_struct<T_y>::value) {
        for (size_t n = 0; n < stan::length(y); ++n)
          operands_and_partials.d_x1[n] *= P;
      }
      if (!is_constant_struct<T_loc>::value) {
        for (size_t n = 0; n < stan::length(mu); ++n)
          operands_and_partials.d_x2[n] *= P;
      }
      if (!is_constant_struct<T_scale>::value) {
        for (size_t n = 0; n < stan::length(sigma); ++n)
          operands_and_partials.d_x3[n] *= P;
      }

      return operands_and_partials.to_var(P, y, mu, sigma);
    }
  }
}
#endif
