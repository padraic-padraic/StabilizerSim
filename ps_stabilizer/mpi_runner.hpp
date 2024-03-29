#ifndef MPIRUNNER_HPP
#define MPIRUNNER_HPP

#include "core.hpp"
#include "rng.hpp"
#include "runner.hpp"

#include "mpi.h"

namespace StabilizerSimulator
{

static std::vector<int> EMPTY_RANKS;

template<class stabilizer_t> class MPIRunner : public Runner<stabilizer_t>
{
public:

  using BaseRunner = Runner<stabilizer_t>;

  MPIRunner(std::vector<int> ranks, MPI_Comm communicator) :
  group(MPI_GROUP_NULL),
  finalized(false),
  is_in_lib(false),
  is_master(false)
  {
    MPI_Comm to_dup = communicator;
    if (to_dup == MPI_COMM_NULL)
    {
      to_dup = MPI_COMM_WORLD;
    }
    if (ranks.size() == 0)
    {
      MPI_Comm_dup(to_dup, &comm);
      MPI_Comm_group(comm, &group);
    }
    else 
    {
      MPI_Group dup_group;
      MPI_Comm_group(to_dup, &dup_group);
      MPI_Group_incl(dup_group, ranks.size(), ranks.data(), &group);
      MPI_Comm_create_group(to_dup, group, 0, &comm);
    }
    if (comm == MPI_COMM_NULL)
    {
      is_in_lib = false;
      my_rank = MPI_UNDEFINED;
      n_procs = MPI_UNDEFINED;
    }
    else
    {
      is_in_lib = true;
      MPI_Comm_size(comm, &n_procs);
      MPI_Comm_rank(comm, &my_rank);
      if (my_rank == 0)
      {
        is_master = true;
      }
    }
  };

  MPIRunner() : MPIRunner(EMPTY_RANKS, MPI_COMM_NULL) {};

  MPIRunner(unsigned n_qubits_, unsigned n_states_) : 
  MPIRunner(EMPTY_RANKS, MPI_COMM_NULL)
  {
    stabilizer_t init(n_qubits_);
    initialize(n_qubits_, n_states_, init);
  };

  void initialize(unsigned n_qubits_, unsigned n_states_, stabilizer_t& initial_state);
  void finalize();

private:
  MPI_Comm comm;
  MPI_Group group;

  bool finalized;
  int n_procs;
  int my_rank;
  bool is_in_lib;
  bool is_master;
  
  unsigned total_states;
  int index_low;  

  int i4_div_rounded(int a, int b);

  void init_metropolis() override;
  void metropolis_step() override;

};

template <class stabilizer_t>
void MPIRunner<stabilizer_t>::initialize(unsigned n_qubits_, unsigned n_states_, stabilizer_t& initial_state)
{
  BaseRunner::n_qubits = n_qubits_;
  total_states = n_states_;
  unsigned seed = std::time(nullptr);
  MPI_Bcast(&seed, 1, MPI_UNSIGNED, 0, comm);
  init_rng(seed, my_rank);
  if (my_rank==0) //Divive up the states
  {
    int task_divider = n_procs;
    int states_chunk = i4_div_rounded(total_states, task_divider);
    BaseRunner::n_states = states_chunk;
    index_low = 0;
    // Setup divider loop
    task_divider= task_divider-1;
    int chi_left = total_states - states_chunk;
    int low_index = states_chunk;
    for (int proc=1; proc<n_procs; proc++)
    {
      unsigned send_states_chunk = i4_div_rounded(chi_left, task_divider);
      MPI_Send(&send_states_chunk, 1, MPI_UNSIGNED, proc, 0, comm);
      MPI_Send(&low_index, 1, MPI_INT, proc, 0, comm);
      chi_left = chi_left - send_states_chunk;
      low_index += send_states_chunk;
      task_divider--;
    }
  }
  else
  {
    MPI_Recv(&BaseRunner::n_states, 1, MPI_UNSIGNED, 0, 0, comm, MPI_STATUS_IGNORE);
    MPI_Recv(&index_low, 1, MPI_INT, 0, 0, comm, MPI_STATUS_IGNORE);
  }
  BaseRunner::states.clear();
  BaseRunner::coefficients.clear();
  BaseRunner::states.reserve(BaseRunner::n_states);
  BaseRunner::coefficients.reserve(BaseRunner::n_states);
  for(unsigned i=0; i<BaseRunner::n_states; i++)
  {
    BaseRunner::states.push_back(initial_state);
    BaseRunner::coefficients.push_back(complex_t(1., 0.));
  }
}

template <class stabilizer_t>
void MPIRunner<stabilizer_t>::finalize()
{
  if (group != MPI_GROUP_NULL)
  {
    MPI_Group_free(&group);
  }
  if (is_in_lib)
  {
    MPI_Comm_free(&comm);    
  }
  finalized = true;
}

template <class stabilizer_t>
int MPIRunner<stabilizer_t>::i4_div_rounded(int a, int b)
{
  // Implementation from https://people.sc.fsu.edu/~jburkardt/cpp_src/task_division/task_division.html
  // Written by  John Burkardt, Department of Scientific Computing, Florida State University
  // Used under the LGPL License
  int a_abs;
  int b_abs;
  static int i4_huge = 2147483647;
  int value;
  if ( a == 0 && b == 0 )
  {
    value = i4_huge;
  }
  else if ( a == 0 )
  {
    value = 0;
  }
  else if ( b == 0 )
  {
    if ( a < 0 )
    {
      value = - i4_huge;
    }
    else
    {
      value = + i4_huge;
    }
  }
  else
  {
    a_abs = abs ( a );
    b_abs = abs ( b );

    value = a_abs / b_abs;
    if ( ( 2 * value + 1 ) * b_abs < 2 * a_abs )
    {
      value = value + 1;
    }
    if ( ( a < 0 && 0 < b ) || ( 0 < a && b < 0 ) )
    {
      value = - value;
    }
  }
  return value;
}

template <class stabilizer_t>
void MPIRunner<stabilizer_t>::init_metropolis()
{
  BaseRunner::accept = 0;
  BaseRunner::last_proposal = 0;
  //Random initial x_string from RngEngine
  if(is_master)
  {
    uint_t max = (1ULL<<BaseRunner::n_qubits) - 1;
    BaseRunner::x_string = (random_uint() & max);
  }
  MPI_Bcast(&BaseRunner::x_string, 1, MPI_UNSIGNED_LONG_LONG, 0, comm);
  double local_real=0., local_imag=0.;
  double real_ampsum=0., imag_ampsum=0.;
  complex_t local(0., 0.);
  for (uint_t i=0; i<BaseRunner::n_states; i++)
  {
    scalar_t amp = BaseRunner::states[i].Amplitude(BaseRunner::x_string);
    if(amp.eps == 1)
    {
      local += (amp.to_complex() * BaseRunner::coefficients[i]);
    }
  }
  local_real = local.real();
  local_imag = local.imag();
  MPI_Reduce(&local_real, &real_ampsum, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
  MPI_Reduce(&local_imag, &imag_ampsum, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
  if(is_master)
  {
    BaseRunner::old_ampsum = complex_t(real_ampsum, imag_ampsum);
  }
}

template <class stabilizer_t>
void MPIRunner<stabilizer_t>::metropolis_step()
{
  MPI_Bcast(&BaseRunner::last_proposal, 1, MPI_UNSIGNED, 0, comm);
  MPI_Bcast(&BaseRunner::accept, 1, MPI_UNSIGNED, 0, comm);
  unsigned proposal;
  if(is_master)
  {
    proposal = random_uint() % BaseRunner::n_qubits; //TODO: Random number getting!
  }
  MPI_Bcast(&proposal, 1, MPI_UNSIGNED, 0, comm);
  if(BaseRunner::accept)
  {
    BaseRunner::x_string ^= (one << BaseRunner::last_proposal);
  }
  double real_part = 0., imag_part = 0.;
  double local_real =0., local_imag = 0.;
  if (BaseRunner::accept == 0)
  {
    for (uint_t i=0; i<BaseRunner::n_states; i++)
    {
      scalar_t amp = BaseRunner::states[i].ProposeFlip(proposal);
      if(amp.eps == 1)
      {
        complex_t local = amp.to_complex() * BaseRunner::coefficients[i];
        local_real += local.real();
        local_imag += local.imag();
      }
    }
  }
  else
  {
    for (uint_t i=0; i<BaseRunner::n_states; i++)
    {
      BaseRunner::states[i].AcceptFlip();
      scalar_t amp = BaseRunner::states[i].ProposeFlip(proposal);
      if(amp.eps == 1)
      {
        complex_t local = amp.to_complex() * BaseRunner::coefficients[i];
        local_real += local.real();
        local_imag += local.imag();
      }
    }
  }
  MPI_Reduce(&local_real, &real_part, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
  MPI_Reduce(&local_imag, &imag_part, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
  if(is_master)
  {
    complex_t ampsum(real_part, imag_part);
    double p_threshold = std::norm(ampsum)/std::norm(BaseRunner::old_ampsum);
    #ifdef  __FAST_MATH__ //isnan doesn't behave well under fastmath, so use absolute tolerance check instead
    if(std::isinf(p_threshold) || std::abs(std::norm(old_ampsum)-0.) < 1e-8)
    #else
    if(std::isinf(p_threshold) || std::isnan(p_threshold))
    #endif
    {
      BaseRunner::accept = 1;
      BaseRunner::old_ampsum = ampsum;
      BaseRunner::last_proposal = proposal; //We try to move away from node with 0 probability.
    }
    else
    {
      double rand = random_double();
      if (rand < p_threshold)
      {
        BaseRunner::accept = 1;
        BaseRunner::old_ampsum = ampsum;
        BaseRunner::last_proposal = proposal;
      }
      else
      {
        BaseRunner::accept = 0;
      }
    }
  }
}

} //end namespace StabilizerSimulator

#endif