#ifndef RUNNER_HPP
#define RUNNER_HPP

#include "core.hpp"
#include "norm_estimation.hpp"
#include "rng.hpp"

#include <ctime>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

// Runner class
// ME/NE functions as in the Qiskit Runner
// Decomposition and probability functions as friend functors (allows state). 
// Run method also overridable

namespace StabilizerSimulator
{

struct DecompositionBuilder
{
  DecompositionBuilder();
  virtual ~DecompositionBuilder() = default;

  virtual void operator()(unsigned rank) = 0;
};

template <class stabilizer_t> class Runner
{
public:
  Runner();
  Runner(unsigned n_qubits_, unsigned n_states_) :
  n_states(n_states_), n_qubits(n_qubits_),
  states(n_states_, stabilizer_t(n_qubits_)),
  coefficients(n_states_, complex_t(1., 0.))
  {};
  ~Runner() = default;

  friend DecompositionBuilder;

  void initialize(unsigned n_qubits_, unsigned n_states_, stabilizer_t& initial_state);

  void build_decomposition(DecompositionBuilder& db_func);
  uint_t monte_carlo_sampler(unsigned mixing_time=7000);
  std::vector<uint_t> monte_carlo_sampler(unsigned n_shots=3000, unsigned mixing_time=7000);
  double norm_estimation(unsigned n_samples=50);
  // double norm_estimation(ProbabilityFunc& p_func, std::vector<pauli_t>& generators, unsigned n_samples=50);

private:
  unsigned n_states;
  unsigned n_qubits;
  std::vector<stabilizer_t> states;
  std::vector<complex_t> coefficients;

  bool accept;
  complex_t old_ampsum;
  uint_t x_string;
  uint_t last_proposal;

  virtual void init_metropolis();
  virtual void metropolis_step();
};

template <class stabilizer_t>
void Runner<stabilizer_t>::initialize(unsigned n_qubits_, unsigned n_states_, stabilizer_t& initial_state)
{
  states.clear();
  coefficients.clear();
  states.reserve(n_states_);
  coefficients.reserve(n_states_);
  n_qubits = n_qubits_;
  complex_t initial_coeff(1., 0.);
  unsigned seed = std::time(nullptr);
  #ifdef _OPENMP
    #pragma parallel
    {
      init_rng(seed, omp_get_thread_num());
    }
  #else
    init_rng(seed, 0);
  #endif
  for(unsigned i=0; i<n_states_; i++)
  {
    states.push_back(initial_state);
    coefficients.push_back(initial_coeff);
  }
}

template <class stabilizer_t>
void Runner<stabilizer_t>::build_decomposition(DecompositionBuilder& db_func)
{
  #pragma omp parallel for
  for(unsigned i=0; i<n_states; i++)
  {
    db_func(i);
  }
}

template <class stabilizer_t>
void Runner<stabilizer_t>::init_metropolis()
{
  accept = 0;
  uint_t max = (1ULL<<n_qubits) - 1;
  x_string = (rand_int() & max);
  last_proposal=0;
  double local_real=0., local_imag=0.;
  #pragma omp parallel for reduction(+:local_real) reduction(+:local_imag)
  for (uint_t i=0; i<n_states; i++)
  {
    scalar_t amp = states[i].Amplitude(x_string);
    if(amp.eps == 1)
    {
      complex_t local = amp.to_complex() * coefficients[i];
      local_real += local.real();
      local_imag += local.imag();      
    }
  }
  old_ampsum = complex_t(local_real, local_imag);
}

template <class stabilizer_t>
void Runner<stabilizer_t>::metropolis_step()
{
  uint_t proposal = rand_int() % n_qubits;
  if(accept)
  {
    x_string ^= (one << last_proposal);
  }
  double real_part = 0.,imag_part =0.;
  if (accept == 0)
  {
    #pragma omp parallel for reduction(+:real_part) reduction(+:imag_part)
    for (uint_t i=0; i<n_states; i++)
    {
      scalar_t amp = states[i].ProposeFlip(proposal);
      if(amp.eps == 1)
      {
        complex_t local = amp.to_complex() * coefficients[i];
        real_part += local.real();
        imag_part += local.imag();
      }
    }
  }
  else
  {
    #pragma omp parallel for reduction(+:real_part) reduction(+:imag_part)
    for (uint_t i=0; i<n_states; i++)
    {
      states[i].AcceptFlip();
      scalar_t amp = states[i].ProposeFlip(proposal);
      if(amp.eps == 1)
      {
        complex_t local = amp.to_complex() * coefficients[i];
        real_part += local.real();
        imag_part += local.imag();
      }
    }
  }
  complex_t ampsum(real_part, imag_part);
  double p_threshold = std::norm(ampsum)/std::norm(old_ampsum);
  #ifdef  __FAST_MATH__ //isnan doesn't behave well under fastmath, so use absolute tolerance check instead
  if(std::isinf(p_threshold) || std::abs(std::norm(old_ampsum)-0.) < 1e-8)
  #else
  if(std::isinf(p_threshold) || std::isnan(p_threshold))
  #endif
  {
    accept = 1;
    old_ampsum = ampsum;
    last_proposal = proposal; //We try to move away from node with 0 probability.
  }
  else
  {
    double rand = rand_double();
    if (rand < p_threshold)
    {
      accept = 1;
      old_ampsum = ampsum;
      last_proposal = proposal;
    }
    else
    {
      accept = 0;
    }
  }
}

template <class stabilizer_t>
uint_t Runner<stabilizer_t>::monte_carlo_sampler(unsigned mixing_time)
{
  init_metropolis();
  for(unsigned i=0; i<mixing_time; i++)
  {
    metropolis_step();
  }
  return x_string;
}

template <class stabilizer_t>
std::vector<uint_t> Runner<stabilizer_t>::monte_carlo_sampler(unsigned n_shots, unsigned mixing_time)
{
  std::vector<uint_t> shots;
  shots.reserve(n_shots);
  shots.push_back(monte_carlo_sampler(mixing_time));
  for(unsigned i=1; i<n_shots; i++)
  {
    metropolis_step();
    shots.push_back(x_string);
  }
  return shots;
}


} //end namespace Stabilizer Simulator

#endif