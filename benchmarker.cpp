#include "ps_stabilizer/chstabilizer.hpp"
#include "ps_stabilizer/dchstabilizer.hpp"
#include "chp/chp.hpp"
#include "GraphSim/graphsim.hpp"

#include <chrono>
#include <cstdlib>
#include <iostream>
#include <map>
#include <fstream>
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
      CHP::x(&state, control);
      break;
  case Gate::Y:
      CHP::y(&state, control);
      break;
  case Gate::Z:
      CHP::z(&state, control);
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

unsigned DEFAULT_REPETITIONS = 1000;

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

enum class Basis {X, Y, Z};

std::map<std::string, Basis> basis_names(
{
  {"X", Basis::X},
  {"Y", Basis::Y},
  {"Z", Basis::Z}
});

//Global Measurement Variables

int CHPRES;
int GSRES;


template<class T> std::chrono::duration<double> benchmark_measure(T &state, unsigned n_qubits, unsigned warmup, Basis b)
{
  //placeholder
  random_circuit(state, ALL_GATES, n_qubits, warmup);
  unsigned qubit = rand() % n_qubits;
  StabilizerSimulator::pauli_t P;
  switch (b)
  {
    case Basis::X:
      P.X ^= (1ULL << qubit);
      break;
    case Basis::Y:
      P.X ^= (1ULL << qubit);
      P.Z ^= (1ULL << qubit);
    case Basis::Z:
      P.Z ^= (1ULL << qubit);
      break;
  }
  auto start = std::chrono::high_resolution_clock::now();
  state.MeasurePauli(P);
  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> diff = end-start;
  return diff;
}

template<> std::chrono::duration<double> benchmark_measure(CHP::QState &state, unsigned n_qubits, unsigned warmup, Basis b)
{
  //placeholder
  random_circuit(state, ALL_GATES, n_qubits, warmup);
  unsigned qubit = rand() % n_qubits;
  StabilizerSimulator::pauli_t P;
  switch (b)
  {
    case Basis::X:
      CHP::hadamard(&state, qubit);
      break;
    case Basis::Y:
      CHP::hadamard(&state, qubit);
      CHP::phase(&state, qubit);
      break;
  }
  auto start = std::chrono::high_resolution_clock::now();
  CHPRES = CHP::measure(&state, qubit, 1);
  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> diff = end-start;
  return diff;
}

template<> std::chrono::duration<double> benchmark_measure(GraphSim::GraphRegister &state, unsigned n_qubits, unsigned warmup, Basis b)
{
  GraphSim::LocCliffOp basis_choice = GraphSim::lco_Z;
  random_circuit(state, ALL_GATES, n_qubits, warmup);
  unsigned qubit = rand() % n_qubits;
  switch(b)
  {
    case Basis::X:
      basis_choice = GraphSim::lco_X;
      break;
    case Basis::Y:
      basis_choice = GraphSim::lco_Y;
      break;
  }
  auto start = std::chrono::high_resolution_clock::now();
  GSRES = state.measure(qubit, basis_choice);
  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> diff = end-start;
  return diff;
}

template<class T> std::chrono::duration<double> benchmark_innerprod(T &state, T &state2, unsigned n_qubits, unsigned warmup)
{
  //placeholder
  random_circuit(state, ALL_GATES, n_qubits, warmup);
  random_circuit(state2, ALL_GATES, n_qubits, warmup);
  auto start = std::chrono::high_resolution_clock::now();
  // state.InnerProduct(state2);
  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> diff = end-start;
  return diff;
}

template<class T> std::chrono::duration<double> benchmark_op(T &state, unsigned n_qubits,
                                                             std::string op_string, unsigned warmup, Basis b)
{
  random_circuit(state, ALL_GATES, n_qubits, warmup);
  auto is_gate = gate_names.find(op_string);
  if(is_gate != gate_names.end())
  {
    return benchmark_gate(state, n_qubits, is_gate->second);
  }
  if(op_string == "measure")
  {
    return benchmark_measure(state, n_qubits, warmup, b);
  }
  throw std::runtime_error("Don't recognise op_string: " + op_string);
}

int main(int argc, char* argv[])
{
  if(argc < 5)
  {
      throw std::runtime_error("Expected at least four arguments: simulator, operation, beta, and an output file name.");
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
  std::ofstream output_data;
  std::string out_name = argv[4];
  output_data.open(out_name);
  output_data << "Num Qubits \t Average Operation Time\n";
  Basis b = Basis::Z;
  if (argc == 6)
  {
    if (op_string != "measure")
    {
      std::cout << "Extra arguments ignored unless this is a measurement" << std::endl;
    }
    else
    {
      std::string basis_choice = argv[5];
      auto basis = basis_names.find(basis_choice);
      if (basis == basis_names.end())
      {
        throw std::runtime_error("Do not recognise measurement basis " + basis_choice);
      }
      b = basis->second;
    }
  }
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
            sum_time += benchmark_op(ch, n, op_string, circuit_warmup, b);
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
            sum_time += benchmark_op(dch, n, op_string, circuit_warmup, b);
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
            sum_time += benchmark_op(chp, n, op_string, circuit_warmup, b);
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
            sum_time += benchmark_op(gsim, n, op_string, circuit_warmup, b);
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
    output_data << n << "\t" << sum_time.count()/DEFAULT_REPETITIONS << "\n";
  }
  output_data.close();
  return 0;
}