#ifndef STAN_MCMC_HMC_NUTS_ADAPT_UNIT_E_NUTS_HPP
#define STAN_MCMC_HMC_NUTS_ADAPT_UNIT_E_NUTS_HPP

#include <stan/mcmc/stepsize_adapter.hpp>
#include <stan/mcmc/hmc/nuts/unit_e_nuts.hpp>

namespace stan {

  namespace mcmc {

    // The No-U-Turn Sampler (NUTS) on a
    // Euclidean manifold with unit metric
    // and adaptive stepsize

    template <typename Model, class BaseRNG>
    class adapt_unit_e_nuts: public unit_e_nuts<Model, BaseRNG>,
                             public stepsize_adapter {
    public:
      adapt_unit_e_nuts(Model &model, BaseRNG& rng,
                        std::ostream* o, std::ostream* e)
        : unit_e_nuts<Model, BaseRNG>(model, rng, o, e) {}

      ~adapt_unit_e_nuts() {}

      sample transition(sample& init_sample) {
        sample s = unit_e_nuts<Model, BaseRNG>::transition(init_sample);

        if (this->adapt_flag_)
          this->stepsize_adaptation_.learn_stepsize(this->nom_epsilon_,
                                                    s.accept_stat());

        return s;
      }

      void disengage_adaptation() {
        base_adapter::disengage_adaptation();
        this->stepsize_adaptation_.complete_adaptation(this->nom_epsilon_);
      }
    };

  }  // mcmc

}  // stan

#endif
