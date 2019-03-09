#include "ps_stabilizer/chstabilizer.hpp"
#include "ps_stabilizer/dchstabilizer.hpp"
#include "chp/chp.hpp"
#include "GraphSim/graphsim.hpp"

#include <chrono>
#include <cstdlib>
#include <iostream>
#include <map>
#include <string>
#include <vector>

enum class Gate {
  X, Y, Z, H, S, Sdag, CZ, CX
};

std::map<std::string, Gate> gate_names( 
{
  {"X", Gate::X},
  {"Y", Gate::Y},
  {"Z", Gate::Z},
  {"H", Gate::H},
  {"S", Gate::S},
  {"Sdag", Gate::Sdag},
  {"CZ", Gate::CZ},
  {"CX", Gate::CX}
});

template<typename T> void apply_gate(T &state, Gate g, unsigned control, unsigned target)
{
  switch (g)
  {
    case Gate::X:
      state.X(control);
      break;
    case Gate::Y:
      state.Y(control);
      break;
    case Gate::Z:
      state.Z(control);
      break;
    case Gate::H:
      state.H(control);
      break;
    case Gate::S:
      state.S(control);
      break;
    case Gate::Sdag:
      state.Sdag(control);
      break;
    case Gate::CZ:
      state.CZ(control, target);
      break;
    case Gate::CX:
      state.CX(control, target);
      break;
    default:
      throw std::logic_error("Wat");
      break;
    }
}

template<> void apply_gate(CHP::QState &state, Gate g, unsigned control, unsigned target)
{
  switch (g)
  {
  case Gate::X: //This might throw the benchmarks off somewhat?
      CHP::hadamard(&state, control);
      CHP::phase(&state, control);
      CHP::phase(&state, control);
      CHP::hadamard(&state, control);
      break;
  case Gate::Y:
      //X
      CHP::hadamard(&state, control);
      CHP::phase(&state, control);
      CHP::phase(&state, control);
      CHP::hadamard(&state, control);
      //Z
      CHP::phase(&state, control);
      CHP::phase(&state, control);
      break;
  case Gate::Z:
      CHP::phase(&state, control);
      CHP::phase(&state, control);
      break;
  case Gate::H:
      CHP::hadamard(&state, control);
      break;
  case Gate::S:
      CHP::phase(&state, control);
      break;
  case Gate::Sdag:
      CHP::phase(&state, control);
      CHP::phase(&state, control);
      CHP::phase(&state, control);
      break;
  case Gate::CZ:
      CHP::hadamard(&state, target);
      CHP::cnot(&state, control, target);
      CHP::hadamard(&state, target);
      break;
  case Gate::CX:
      CHP::cnot(&state, control, target);
      break;
    default:
      throw std::logic_error("Wat");
      break;
  }
}

template<> void apply_gate(GraphSim::GraphRegister &state, Gate g, unsigned control, unsigned target)
{
  switch (g)
  {
    case Gate::X: //This might throw the benchmarks off somewhat?
      state.local_op(control, GraphSim::lco_X);
      break;
    case Gate::Y:
      state.local_op(control, GraphSim::lco_Y);
      break;
    case Gate::Z:
      state.local_op(control, GraphSim::lco_X);
      break;
    case Gate::H:
      state.hadamard(control);
      break;
    case Gate::S:
      state.local_op(control, GraphSim::lco_S);
      break;
    case Gate::Sdag:
      state.local_op(control, GraphSim::lco_Sh);
      break;
    case Gate::CZ:
      state.cphase(control, target);
      break;
    case Gate::CX:
      state.cnot(control, target);
      break;
    default:
      throw std::logic_error("Wat");
      break;
  }
}

template<class T> void random_circuit(T &state, std::vector<Gate> &gates, unsigned n_qubits, unsigned n_gates=25)
{
  for(unsigned i=0; i<n_gates; i++)
  {
    unsigned control = (rand()%n_qubits);
    unsigned target = (rand()%n_qubits);
    if(target == control)
    {
        target = (target +1)%n_qubits;
    }
    unsigned gate_index = (rand()%gates.size());
    Gate g = gates[gate_index];
    apply_gate(state, g, control, target);
  }
}

unsigned DEFAULT_REPETITIONS = 1000000;

enum class Representation
{
  CH,
  DCH,
  CHP,
  GRAPHSIM
};

std::map<std::string, Representation> rep_names(
{
  {"ch", Representation::CH},
  {"dch", Representation::DCH},
  {"chp", Representation::CHP},
  {"graphsim", Representation::GRAPHSIM}
});

std::vector<Gate> ALL_GATES = {Gate::X, Gate::Y, Gate::Z, Gate::H, Gate::S, Gate::Sdag, Gate::CZ, Gate::CX};

template<class T> std::chrono::duration<double> benchmark_gate(T &state, unsigned n_qubits, Gate g)
{
  unsigned control = (rand()%n_qubits);
  unsigned target = (rand()%n_qubits);
  if(target == control)
  {
    target = (target +1)%n_qubits;
  }
  auto start = std::chrono::high_resolution_clock::now();
  apply_gate(state, g, control, target);
  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> diff = end-start;
  return diff;
}

template<class T> std::chrono::duration<double> benchmark_measure(T &state, unsigned n_qubits)
{
  //placeholder
  return std::chrono::duration<double>(0);
}

template<class T> std::chrono::duration<double> benchmark_innerprod(T &state, T &state2, unsigned n_qubits, unsigned warmup)
{
  //placeholder
  random_circuit(state, ALL_GATES, n_qubits, warmup);
  random_circuit(state, ALL_GATES, n_qubits, warmup);
  auto start = std::chrono::high_resolution_clock::now();
  // state.InnerProduct(state2);
  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> diff = end-start;
  return std::chrono::duration<double>(0);
}

template<class T> std::chrono::duration<double> benchmark_op(T &state, unsigned n_qubits,
                                                             std::string op_string, unsigned warmup)
{
  random_circuit(state, ALL_GATES, n_qubits, warmup);
  auto is_gate = gate_names.find(op_string);
  if(is_gate != gate_names.end())
  {
    return benchmark_gate(state, n_qubits, is_gate->second);
  }
  if(op_string == "measure")
  {
    return benchmark_measure(state, n_qubits);
  }
  throw std::runtime_error("Don't recognise op_string: " + op_string);
}

int main(int argc, char* argv[])
{
  if(argc < 3)
  {
      throw std::runtime_error("Expected at least three arguments: simulator, operation, and beta.");
  }
  std::string sim_name = argv[1];
  auto sim = rep_names.find(sim_name);
  if(sim == rep_names.end())
  {
      throw std::runtime_error("Do not recognise simulator " + sim_name);
  }
  std::string op_string = argv[2];
  double beta = std::stod(argv[3]);
  std::cout << "Running with Simulator: " << sim_name << " and beta: " << beta << std::endl;
  for(unsigned n=5; n<62; n++)
  {
    std::cout << n << " qubits run." << std::endl;
    unsigned circuit_warmup = std::lrint(beta*n*std::log2(n));
    std::chrono::duration<double> sum_time(0);
    for(unsigned i=0; i<DEFAULT_REPETITIONS; i++)
    {
      switch(sim->second)
      {
        case Representation::CH:
        {
          StabilizerSimulator::CHState ch(n);
          if(op_string != "innerprod")
          {
            sum_time += benchmark_op(ch, n, op_string, circuit_warmup);
          }
          else
          {
            StabilizerSimulator::CHState ch2(n);
            sum_time += benchmark_innerprod(ch, ch2, n, circuit_warmup);
          }
          break;
        }
        case Representation::DCH:
        {
          StabilizerSimulator::DCHState dch(n);
          if(op_string != "innerprod")
          {
            sum_time += benchmark_op(dch, n, op_string, circuit_warmup);
          }
          else
          {
            StabilizerSimulator::DCHState dch2(n);
            sum_time += benchmark_innerprod(dch, dch2, n, circuit_warmup);
          }
          break;
        }
        case Representation::CHP:
        {
          CHP::QState chp;
          CHP::initstae_(&chp, n, NULL);
          if(op_string != "innerprod")
          {
            sum_time += benchmark_op(chp, n, op_string, circuit_warmup);
          }
          else
          {
            CHP::QState chp2;
            CHP::initstae_(&chp2, n, NULL);
            sum_time += benchmark_innerprod(chp, chp2, n, circuit_warmup);
          }
          break;
        }
        case Representation::GRAPHSIM:
        {
          GraphSim::GraphRegister gsim(n);
          if(op_string != "innerprod")
          {
            sum_time += benchmark_op(gsim, n, op_string, circuit_warmup);
          }
          else
          {
            throw std::runtime_error("Graphsim doesn't support the Inner Product operation.");
          }
          break;
        }
      }
    }
    std::cout << "Cumulative time: " << sum_time.count() << "Average: " << sum_time.count()/DEFAULT_REPETITIONS << std::endl; 
  }
}