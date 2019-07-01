#include "core.hpp"

#include <complex>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace StabilizerSimulator
{

template<class T> double NormEstimate(std::vector<T>& states, const std::vector< std::complex<double> >& phases, 
                  const std::vector<uint_t>& Samples_d1, const std::vector<uint_t> &Samples_d2, 
                  const std::vector< std::vector<uint_t> >& Samples)
{
  int n_threads = 1;
  #ifdef _OPENMP
  n_threads = omp_get_max_threads();
  #endif
  return NormEstimate(states, phases, Samples_d1, Samples_d2, Samples, n_threads);
}

template<class T> double NormEstimate(std::vector<T>& states, const std::vector< std::complex<double> >& phases, 
                    const std::vector<uint_t>& Samples_d1, const std::vector<uint_t> &Samples_d2, 
                    const std::vector< std::vector<uint_t> >& Samples, int n_threads)
{
    #ifdef _OPENMP
    if(n_threads == -1)
    {
      n_threads = omp_get_max_threads();
    }
    #endif
    // Norm estimate for a state |psi> = \sum_{i} c_{i}|phi_{i}>
    double xi=0;
    unsigned L = Samples_d1.size();
    // std::vector<double> data = (L,0.);
    for (size_t i=0; i<L; i++)
    {
      double re_eta =0., im_eta = 0.;
      #pragma omp parallel for if(n_threads > 1) reduction(+:re_eta) reduction(+:im_eta) num_threads(n_threads)
      for (uint_t j=0; j<states.size(); j++)
      {
          if(states[j].ScalarPart().eps != 0)
          {
              scalar_t amp = states[j].InnerProduct(Samples_d1[i], Samples_d2[i], Samples[i]);
              if (amp.eps != 0)
              {
                  if (amp.e % 2)
                  {
                    amp.p--;
                  }
                  double mag = pow(2, amp.p/(double) 2);
                  std::complex<double> phase(RE_PHASE[amp.e], IM_PHASE[amp.e]);
                  phase *= conj(phases[j]);
                  re_eta += (mag * real(phase));
                  im_eta += (mag * imag(phase));
              }
          }
      }
      xi += (pow(re_eta,2) + pow(im_eta, 2));
    }
    return std::pow(2, states[0].NQubits()) * (xi/L);
}

}
