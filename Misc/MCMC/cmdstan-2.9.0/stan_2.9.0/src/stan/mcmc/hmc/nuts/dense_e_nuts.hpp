#ifndef STAN_MCMC_HMC_NUTS_DENSE_E_NUTS_HPP
#define STAN_MCMC_HMC_NUTS_DENSE_E_NUTS_HPP

#include <stan/mcmc/hmc/nuts/base_nuts.hpp>
#include <stan/mcmc/hmc/hamiltonians/dense_e_point.hpp>
#include <stan/mcmc/hmc/hamiltonians/dense_e_metric.hpp>
#include <stan/mcmc/hmc/integrators/expl_leapfrog.hpp>

namespace stan {
  namespace mcmc {

    // The No-U-Turn Sampler (NUTS) on a
    // Euclidean manifold with dense metric
    template <typename Model, class BaseRNG>
    class dense_e_nuts : public base_nuts<Model, dense_e_metric,
                                          expl_leapfrog, BaseRNG> {
    public:
      dense_e_nuts(Model &model, BaseRNG& rng, std::ostream* o,
                   std::ostream* e)
        : base_nuts<Model, dense_e_metric, expl_leapfrog,
                    BaseRNG>(model, rng, o, e) {
        this->name_ = "NUTS with a dense Euclidean metric";
      }

      // Note that the points don't need to be swapped
      // here since start.mInv = finish.mInv
      bool compute_criterion(ps_point& start,
                             dense_e_point& finish,
                             Eigen::VectorXd& rho) {
        return finish.p.transpose() * finish.mInv * (rho - finish.p) > 0
               && start.p.transpose() * finish.mInv * (rho - start.p)  > 0;
      }
    };

  }  // mcmc
}  // stan
#endif
