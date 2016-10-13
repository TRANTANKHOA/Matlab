#ifndef STAN_MCMC_HMC_HAMILTONIANS_DIAG_E_METRIC_HPP
#define STAN_MCMC_HMC_HAMILTONIANS_DIAG_E_METRIC_HPP

#include <boost/random/variate_generator.hpp>
#include <boost/random/normal_distribution.hpp>

#include <stan/mcmc/hmc/hamiltonians/base_hamiltonian.hpp>
#include <stan/mcmc/hmc/hamiltonians/diag_e_point.hpp>

namespace stan {

  namespace mcmc {

    // Euclidean manifold with diagonal metric
    template <typename Model, typename BaseRNG>
    class diag_e_metric: public base_hamiltonian<Model, diag_e_point, BaseRNG> {
    public:
      diag_e_metric(Model& model, std::ostream* e)
        : base_hamiltonian<Model, diag_e_point, BaseRNG>(model, e) {}

      ~diag_e_metric() {}

      double T(diag_e_point& z) {
        return 0.5 * z.p.dot( z.mInv.cwiseProduct(z.p) );
      }

      double tau(diag_e_point& z) {
        return T(z);
      }

      double phi(diag_e_point& z) {
        return this->V(z);
      }

      const Eigen::VectorXd dtau_dq(diag_e_point& z) {
        return Eigen::VectorXd::Zero(this->model_.num_params_r());
      }

      const Eigen::VectorXd dtau_dp(diag_e_point& z) {
        return z.mInv.cwiseProduct(z.p);
      }

      const Eigen::VectorXd dphi_dq(diag_e_point& z) {
        return z.g;
      }

      void sample_p(diag_e_point& z, BaseRNG& rng) {
        boost::variate_generator<BaseRNG&, boost::normal_distribution<> >
          rand_diag_gaus(rng, boost::normal_distribution<>());

        for (int i = 0; i < z.p.size(); ++i)
          z.p(i) = rand_diag_gaus() / sqrt(z.mInv(i));
      }
    };

  }  // mcmc

}  // stan

#endif
