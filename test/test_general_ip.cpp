#define CATCH_CONFIG_RUNNER
#include "test/lib/catch.hpp"

#include "ps_stabilizer/core.hpp"
#include "ps_stabilizer/chstabilizer.hpp"
#include "ps_stabilizer/dchstabilizer.hpp"

#include <algorithm>
#include <cstdlib>
#include <complex>
#include <iostream>
#include <map>
#include <string>
#include <time.h>
#include <vector>

using namespace StabilizerSimulator;

using complex_t = std::complex<double>;

//--------------------------------//
// Helper methods                 //
//--------------------------------//

const double precision=1e-8;

bool complex_close(complex_t &a, complex_t &b)
{
    return( (std::abs(a.real()-b.real())<precision) & (std::abs(a.imag()-b.imag())<precision) );
}

bool check_states(CHState &ch, DCHState &dch)
{
    unsigned n = ch.NQubits();
    uint_t dim = one << n;
    bool result = true;
    for(uint_t i=0; i<dim; i++)
    {
        scalar_t amp;
        complex_t ch_amp, dch_amp;
        amp = ch.Amplitude(i);
        ch_amp = amp.to_complex();
        amp = dch.Amplitude(i);
        dch_amp = amp.to_complex();
        if(!complex_close(ch_amp, dch_amp))
        {
            std::cout << ch_amp << " != " << dch_amp << " for string ";
            Print(i, n);
            std::cout << std::endl;
            result =  false;
        }
    }
    return result;
}

template <class T> std::vector<complex_t> svector(T& state)
{
    uint_t dim = (one << state.n);

}

bool check_states(CHState &ch, CHState &ch2)
{
    unsigned n = ch.NQubits();
    uint_t dim = one << n;
    bool result = true;
    for(uint_t i=0; i<dim; i++)
    {
        scalar_t amp;
        complex_t ch_amp, ch2_amp;
        amp = ch.Amplitude(i);
        ch_amp = amp.to_complex();
        amp = ch2.Amplitude(i);
        ch2_amp = amp.to_complex();
        if(!complex_close(ch_amp, ch2_amp))
        {
            std::cout << ch_amp << " != " << ch2_amp << " for string ";
            Print(i, n);
            std::cout << std::endl;
            result =  false;
        }
    }
    return result;
}

//--------------------------------//
// Generate random circuits       //
//--------------------------------//

enum class Gate {
    X, Y, Z, H, S, Sdag, CZ, CX
};

std::map<Gate, std::string> gate_names( 
{
    {Gate::X, "X"},
    {Gate::Y, "Y"},
    {Gate::Z, "Z"},
    {Gate::H, "H"},
    {Gate::S, "S"},
    {Gate::Sdag, "Sdag"},
    {Gate::CZ, "CZ"},
    {Gate::CX, "CX"}
});

std::map<std::string, Gate> gate_enums( 
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

std::map<std::string, Gate> conjugate_enums( 
{
  {"X", Gate::X},
  {"Y", Gate::Y},
  {"Z", Gate::Z},
  {"H", Gate::H},
  {"S", Gate::Sdag},
  {"Sdag", Gate::S},
  {"CZ", Gate::CZ},
  {"CX", Gate::CX}
});

template<class T> void apply_gate(T &state, Gate g, unsigned control, unsigned target, bool echo=false)
{
    if(echo)
    {
        std::cout << "Applying: " << gate_names[g] << "(" << control << "," << target << ")" << std::endl;
    }
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

template<class T> std::string random_circuit(T &state, std::vector<Gate> &gates, unsigned n_gates=25)
{
    unsigned n_qubits = state.NQubits();
    std::string gate_seq = "";
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
        auto gate_str = gate_names[g];
        gate_str += "(" + std::to_string(control);
        if (g==Gate::CZ || g==Gate::CX)
        {
            gate_str+= "," + std::to_string(target);
        }
        gate_str+=")";
        gate_seq += gate_str;
    }
    return gate_seq;
}

void random_circuit_verify(CHState &state1,  DCHState &state2, std::vector<Gate> &gates, unsigned n_gates=25)
{
    REQUIRE(state1.NQubits() == state2.NQubits());
    unsigned n_qubits = state1.NQubits();
    std::string gate_seq = "";
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
        apply_gate(state1, g, control, target);
        apply_gate(state2, g, control, target);
        auto gate_str = gate_names[g];
        gate_str += "(" + std::to_string(control);
        if (g==Gate::CZ || g==Gate::CX)
        {
            gate_str+= "," + std::to_string(target);
        }
        gate_str+=")_";
        gate_seq += gate_str;
        // std::vector<uint_t> M = state2.MMatrix();
        // for(unsigned i=0; i<n_qubits; i++)
        // {
        //     if(!!(M[i] & (one << i)))
        //     {
        //         std::cout << "Found non-zero diagonal in M after doing " << gate_str << std::endl;
        //         REQUIRE(false);
        //     }
        // }
        INFO("Gates: " << gate_seq);
        REQUIRE(check_states(state1, state2));
    }
}

//--------------------------------//
// Circuit types                  //
//--------------------------------//

std::vector<Gate> phase_gates = {Gate::S, Gate::Sdag, Gate::Z, Gate::CZ};
std::vector<Gate> pauli_gates = {Gate::X, Gate::Y, Gate::Z};
std::vector<Gate> deterministic_gates = {Gate::Z, Gate::CZ, Gate::CX, Gate::S, Gate::Sdag, Gate::X, Gate::Y};
std::vector<Gate> all_gates = {Gate::X, Gate::Y, Gate::Z, Gate::S, Gate::Sdag, Gate::CX, Gate::CZ, Gate::H};
std::vector<Gate> c_gates = {Gate::Z, Gate::S, Gate::Sdag, Gate::CZ, Gate::CX};
//--------------------------------//
// Tests                          //
//--------------------------------//

TEST_CASE("Tableu combining")
{
    unsigned n_qubits = 10;
    unsigned n_gates = 2000;
    CHState ch1(n_qubits);
    CHState ch2(n_qubits);
    random_circuit(ch1, all_gates);
    std::string circuit_2 = random_circuit(ch2, c_gates, n_gates);
    std::vector<uint_t> left_G = ch2.GMatrix();
    std::vector<uint_t> left_F = ch2.FMatrix();
    std::vector<uint_t> left_M = ch2.MMatrix();
    uint_t left_g1 = ch2.Gamma1();
    uint_t left_g2 = ch2.Gamma2();
    std::vector<uint_t> right_G = ch1.GMatrix();
    std::vector<uint_t> right_F = ch1.FMatrix();
    std::vector<uint_t> right_M = ch1.MMatrix();    
    std::vector<uint_t> right_GT(n_qubits, zer);
    std::vector<uint_t> right_FT(n_qubits, zer);
    std::vector<uint_t> right_MT(n_qubits, zer);
    for(unsigned i=0; i<n_qubits; i++)
    {
        uint_t shift = (one << i);
        for(unsigned j=0; j<n_qubits; j++)
        {
            if(!!(right_G[j] & shift))
            {
                right_GT[i] ^= (one << j);
            }
            if(!!(right_F[j] & shift))
            {
                right_FT[i] ^= (one << j);
            }
            if(!!(right_M[j] & shift))
            {
                right_MT[i] ^= (one << j);
            }
        }
    }
    uint_t right_g1 = ch1.Gamma1();
    uint_t right_g2 = ch1.Gamma2();
    std::vector<uint_t> new_GT(n_qubits, zer);
    std::vector<uint_t> new_FT(n_qubits, zer);
    std::vector<uint_t> new_MT(n_qubits, zer);
    uint_t new_g1 = zer;
    uint_t new_g2 = zer;
    new_g2 = left_g2;
    new_g1 = left_g1;
    for(unsigned i=0; i<n_qubits; i++)
    {
        uint_t shift = (one << i);
        pauli_t X_out;
        for(unsigned j=0; j<n_qubits; j++)
        {
            //Udpdate the Z stabilizers
            if(!!(left_G[j] & shift))
            {
                new_GT[i] ^= right_GT[j];
            }
            if(!!(left_F[j] & shift))
            {
                pauli_t p;
                p.e += 1*((right_g1>>j) & one);
                p.e += 2*((right_g2>>j) & one);
                p.X = right_FT[j];
                p.Z = right_MT[j];
                X_out *= p;
            }
        }
        new_FT[i] = X_out.X;
        new_MT[i] = X_out.Z;
        for(unsigned j=0; j<n_qubits; j++)
        {
            if(!!(left_M[j] & shift))
            {
                new_MT[i] ^= right_GT[j];
            }
        }
        if(X_out.e & 2)
        {
            new_g2 ^= shift;
        }
        if(X_out.e & 1)
        {
            new_g2 ^= (new_g1&shift);
            new_g1 ^= shift;
        }
    }
    std::vector<uint_t> new_G(n_qubits, zer);
    std::vector<uint_t> new_M(n_qubits, zer);
    std::vector<uint_t> new_F(n_qubits, zer);
    for(unsigned i=0; i<n_qubits; i++)
    {
        for(unsigned j=0; j<n_qubits; j++)
        {
            if(!!(new_GT[i] & (one << j)))
            {
                new_G[j] ^= (one << i);
            }
            if(!!(new_FT[i] & (one << j)))
            {
                new_F[j] ^= (one << i);
            }
            if(!!(new_MT[i] & (one << j)))
            {
                new_M[j] ^= (one << i);
            }            
        }
    }
    ch2.SetG(new_G);
    ch2.SetF(new_F);
    ch2.SetM(new_M);
    ch2.SetGamma1(new_g1);
    ch2.SetGamma2(new_g2);
    //Apply the sequence to right state
    for(unsigned i=0; i<n_gates; i++)
    {
        auto pos = circuit_2.find_first_of("(");
        std::string gstr = circuit_2.substr(0, pos);
        auto found = gate_enums.find(gstr);
        if(found == gate_enums.end())
        {
            throw std::runtime_error("Couldn't find gate " + gstr);
        }
        unsigned control = 0, target = 0;
        if(found->second == Gate::CX || found->second == Gate::CZ)
        {
            auto pos2 = circuit_2.find_first_of(",");
            control = std::stoul(circuit_2.substr(pos2-1, 1));
            target = std::stoul(circuit_2.substr(pos2+1, 1));
        }
        else
        {
            control = std::stoul(circuit_2.substr(pos+1, 1));
        }
        auto end_pos = circuit_2.find(")");
        circuit_2.erase(0, end_pos+1);
        apply_gate(ch1, found->second, control, target);
    }
    right_G = ch1.GMatrix();
    right_F = ch1.FMatrix();
    right_M = ch1.MMatrix();
    for(unsigned i=0; i<n_qubits; i++)
    {
        CHECK(right_G[i] == new_G[i]);
        CHECK(right_F[i] == new_F[i]);
        CHECK(right_M[i] == new_M[i]);
    }
    CHECK(new_g1 == ch1.Gamma1());
    CHECK(new_g2 == ch1.Gamma2());
}

TEST_CASE("Conjugate Tableu combining")
{
    unsigned n_qubits = 10;
    unsigned n_gates = 20;
    CHState ch2(n_qubits);
    uint_t x_string = zer;
    for(unsigned i=0; i<n_qubits; i++)
    {
        if(rand()%2)
        {
            x_string ^= (one << i);
        }
    }
    ch2.CompBasisVector(x_string);
    std::string circuit_2 = random_circuit(ch2, c_gates, n_gates);
    std::cout << "===== SEQUENCE =====" << std::endl;
    std::cout << circuit_2 << std::endl;
    CHState ch1(ch2);
    std::vector<uint_t> left_G = ch2.GMatrix();
    std::vector<uint_t> left_F = ch2.FMatrix();
    std::vector<uint_t> left_M = ch2.MMatrix();
    std::vector<uint_t> left_FT(n_qubits, zer);
    std::vector<uint_t> left_MT(n_qubits, zer);
    uint_t left_g1 = ch2.Gamma1();
    uint_t left_g2 = ch2.Gamma2();
    std::vector<uint_t> right_G = ch1.GMatrix();
    std::vector<uint_t> right_F = ch1.FMatrix();
    std::vector<uint_t> right_M = ch1.MMatrix();    
    std::vector<uint_t> right_GT(n_qubits, zer);
    std::vector<uint_t> right_FT(n_qubits, zer);
    std::vector<uint_t> right_MT(n_qubits, zer);
    for(unsigned i=0; i<n_qubits; i++)
    {
        uint_t shift = (one << i);
        for(unsigned j=0; j<n_qubits; j++)
        {
            if(!!(right_G[j] & shift))
            {
                right_GT[i] ^= (one << j);
            }
            if(!!(right_F[j] & shift))
            {
                right_FT[i] ^= (one << j);
            }
            if(!!(right_M[j] & shift))
            {
                right_MT[i] ^= (one << j);
            }
            if(!!(left_F[j] & shift))
            {
                left_FT[i] ^= (one << j);
            }
            if(!!(left_M[j] & shift))
            {
                left_MT[i] ^= (one << j);
            }
        }
    }
    uint_t right_g1 = ch1.Gamma1();
    uint_t right_g2 = ch1.Gamma2();
    std::vector<uint_t> new_GT(n_qubits, zer);
    std::vector<uint_t> new_FT(n_qubits, zer);
    std::vector<uint_t> new_MT(n_qubits, zer);
    uint_t new_g1 = zer;
    uint_t new_g2 = zer;
    new_g2 = left_g2 ^ left_g1;
    new_g1 = left_g1;
    for(unsigned i=0; i<n_qubits; i++)
    {
        new_g2 ^= hamming_parity(left_FT[i]&left_MT[i])*(one << i);
    }
    for(unsigned i=0; i<n_qubits; i++)
    {
        uint_t shift = (one << i);
        pauli_t X_out;
        for(unsigned j=0; j<n_qubits; j++)
        {
            //Udpdate the Z stabilizers
            if(!!(left_G[j] & shift))
            {
                new_GT[i] ^= right_GT[j];
            }
            if(!!(left_F[j] & shift))
            {
                pauli_t p;
                p.e += 1*((right_g1>>j) & one);
                p.e += 2*((right_g2>>j) & one);
                p.X = right_FT[j];
                p.Z = right_MT[j];
                X_out *= p;
            }
        }
        new_FT[i] = X_out.X;
        new_MT[i] = X_out.Z;
        for(unsigned j=0; j<n_qubits; j++)
        {
            if(!!(left_M[j] & shift))
            {
                new_MT[i] ^= right_GT[j];
            }
        }
        if(X_out.e & 2)
        {
            new_g2 ^= shift;
        }
        if(X_out.e & 1)
        {
            new_g2 ^= (new_g1&shift);
            new_g1 ^= shift;
        }
    }
    std::vector<uint_t> new_G(n_qubits, zer);
    std::vector<uint_t> new_M(n_qubits, zer);
    std::vector<uint_t> new_F(n_qubits, zer);
    for(unsigned i=0; i<n_qubits; i++)
    {
        for(unsigned j=0; j<n_qubits; j++)
        {
            if(!!(new_GT[i] & (one << j)))
            {
                new_G[j] ^= (one << i);
            }
            if(!!(new_FT[i] & (one << j)))
            {
                new_F[j] ^= (one << i);
            }
            if(!!(new_MT[i] & (one << j)))
            {
                new_M[j] ^= (one << i);
            }            
        }
    }
    ch2.SetG(new_G);
    ch2.SetF(new_F);
    ch2.SetM(new_M);
    ch2.SetGamma1(new_g1);
    ch2.SetGamma2(new_g2);
    //Apply the sequence to right state
    for(unsigned i=0; i<n_gates; i++)
    {
        auto pos = circuit_2.find_last_of("(");
        auto start_pos = circuit_2.find_last_of(")", circuit_2.length()-2);
        std::string gstr = circuit_2.substr(start_pos+1, pos-start_pos-1);
        auto found = conjugate_enums.find(gstr);
        if(found == conjugate_enums.end())
        {
            throw std::runtime_error("Couldn't find gate " + gstr);
        }
        unsigned control = 0, target = 0;
        if(found->second == Gate::CX || found->second == Gate::CZ)
        {
            auto pos2 = circuit_2.find_last_of(",");
            control = std::stoul(circuit_2.substr(pos2-1, 1));
            target = std::stoul(circuit_2.substr(pos2+1, 1));
        }
        else
        {
            control = std::stoul(circuit_2.substr(pos+1, 1));
        }
        auto end_pos = circuit_2.find_last_of(")");
        circuit_2.erase(start_pos+1);
        apply_gate(ch2, found->second, control, target);
    }
    CHECK(check_states(ch1, ch2));
}

int main( int argc, char* argv[] ) {
  time_t t;
  unsigned seed = (unsigned) time(&t);
  std::cout << "Initialized with seed: " << seed << std::endl;
  srand(seed);

  int result = Catch::Session().run( argc, argv );

  return result;
}
