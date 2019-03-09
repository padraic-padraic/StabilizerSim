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

template<class T> void apply_gate(T &state, Gate g, unsigned control, unsigned target)
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
        gate_str+=")";
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

//--------------------------------//
// Tests                          //
//--------------------------------//

void DCHState::test_commute_pauli()
{
    uint_t x_string = zer, z_string = zer;
    for(unsigned i=0; i<n; i++)
    {
        if(rand()%2)
        {
            x_string ^= (one << i);
        }
        if(rand()%2)
        {
            z_string ^= (one << i);
        }
    }
    pauli_t z_only;
    z_only.Z = z_string;
    CommutePauli(z_only);
    REQUIRE(z_only.X == zer);
}

TEST_CASE("Test initialisation and Amplitude")
{
    unsigned n_qubits = 5;
    CHState ch(n_qubits);
    DCHState dch(n_qubits);
    REQUIRE(check_states(ch, dch));
}

TEST_CASE("Test Computational state initialisation")
{
    unsigned n_qubits = 10;
    SECTION("With X Gates")
    {
        CHState ch(n_qubits);
        DCHState dch(n_qubits);
        for(unsigned i=0; i<10; i++)
        {
            if (rand()%2)
            {
                ch.X(i);
                dch.X(i);
            }
        }
        REQUIRE(check_states(ch, dch));
    }
    SECTION("With CompBasisVector")
    {
        CHState ch(n_qubits);
        DCHState dch(n_qubits);
        uint_t x_init = zer;
        for(unsigned i=0; i<10; i++)
        {
            if (rand()%2)
            {
                x_init ^= one << i;
            }
        }
        ch.CompBasisVector(x_init);
        dch.CompBasisVector(x_init);
        REQUIRE(check_states(ch, dch));
    }
}

TEST_CASE("Test Pauli Gates")
{

    SECTION("Test Z phases")
    {
        unsigned n_qubits = 10;
        CHState ch(n_qubits);
        DCHState dch(n_qubits);
        uint_t x = zer;
        for(unsigned i=0; i<n_qubits; i++)
        {
            if (rand()%2)
            {
                ch.X(i);
                dch.X(i);
                x ^= (one << i);
            }
        }
        for(uint_t i=0; i<n_qubits; i++)
        {
            if((x>>i)&one)
            {
                ch.Z(i);
                dch.Z(i);
                REQUIRE(check_states(ch, dch));
            }
        }
    }

    SECTION("Test Pauli Y")
    {
        unsigned n_qubits = 10;
        CHState ch(n_qubits);
        DCHState dch(n_qubits);
        for(unsigned i=0; i<20; i++)
        {
            unsigned target = (rand() % 10);
            INFO("Target is " << target);
            ch.Y(target);
            dch.Y(target);
            REQUIRE(check_states(ch, dch));
        }
    }

    SECTION("Random Pauli Circuit")
    {
        unsigned n_qubits = 10;
        CHState ch(n_qubits);
        DCHState dch(n_qubits);
        uint_t x_init = zer;
        for(unsigned i=0; i<n_qubits; i++)
        {
            if(rand() %2)
            {
                x_init ^= (one << i);
            }
        }
        ch.CompBasisVector(x_init);
        dch.CompBasisVector(x_init);
        random_circuit_verify(ch, dch, pauli_gates);
    }
}

TEST_CASE("Test S Gates")
{
    unsigned n_qubits = 10;
    CHState ch(n_qubits);
    DCHState dch(n_qubits);
    uint_t x = zer;
    for(unsigned i=0; i<n_qubits; i++)
    {
        if (rand()%2)
        {
            ch.X(i);
            dch.X(i);
            x ^= (one << i);
        }
    }
    for(uint_t i=0; i<n_qubits; i++)
    {
        if((x>>i)&one)
        {
            ch.S(i);
            dch.S(i);
            REQUIRE(check_states(ch, dch));
        }
    }
    for(uint_t i=0; i<n_qubits; i++)
    {
        if((x>>i)&one)
        {
            ch.Sdag(i);
            dch.Sdag(i);
            REQUIRE(check_states(ch, dch));
        }
    }
}

TEST_CASE("Test CZ gate")
{
    unsigned n_qubits = 10;
    CHState ch(n_qubits);
    DCHState dch(n_qubits);
    unsigned target = (rand()%8) + 1;
    ch.X(target);
    dch.X(target);
    ch.CZ(target, 0);
    dch.CZ(0, target);
    std::vector<uint_t> M = dch.MMatrix();
    REQUIRE(!!(M[target] & (one << 0)));
    REQUIRE(!!(M[0] & (one << target)));
    REQUIRE(check_states(ch, dch));
    ch.X(target);
    dch.X(target);
    REQUIRE(check_states(ch, dch));
}

TEST_CASE("Test Commute Pauli")
{
    unsigned n_qubits = 10;
    DCHState dch(n_qubits);
    std::string gate_seq = random_circuit(dch, phase_gates, 10);
    INFO("GATES WERE: " << gate_seq);
    dch.test_commute_pauli();
}

TEST_CASE("Random Phase Circuit")
{
    unsigned n_qubits = 10;
    CHState ch(n_qubits);
    DCHState dch(n_qubits);
    uint_t x_init = zer;
    for(unsigned i=0; i<n_qubits; i++)
    {
        if(rand()%2)
        {
            x_init ^= (one << i);
        }
    }
    ch.CompBasisVector(x_init);
    dch.CompBasisVector(x_init);
    random_circuit_verify(ch, dch, phase_gates);
}

TEST_CASE("Test CX gate")
{
    unsigned n_qubits = 10;
    SECTION("Check computational strings")
    {
        CHState ch(n_qubits);
        DCHState dch(n_qubits);
        unsigned control = (rand()%8);
        INFO("Picked Control: " << control);
        //Check non-trivial CX gate
        ch.X(control);
        dch.X(control);
        ch.CX(control, 9);
        dch.CX(control, 9);
        REQUIRE(check_states(ch, dch));
        //Check we can reset CX action w/ X
        dch.X(9);
        ch.X(9);
        REQUIRE(check_states(ch, dch));
        // Check CX is self inverse
        ch.X(9);
        dch.X(9);
        ch.CX(control, 9);
        dch.CX(control, 9);
        ch.CX(control, 9);
        dch.CX(control, 9);
        REQUIRE(check_states(ch, dch));
        // Check trivial CX
        ch.CX(8, control);
        dch.CX(8, control);
        REQUIRE(check_states(ch, dch));
    }
    SECTION("Check phases")
    {
        CHState ch(n_qubits);
        DCHState dch(n_qubits);
        unsigned target = (rand()%7)+2;
        INFO("Target was " << target);
        ch.S(target);
        dch.S(target);
        ch.CZ(1, target);
        dch.CZ(1, target);
        ch.CX(1, target);
        dch.CX(1, target);
        dch.CX(0, target);
        REQUIRE(check_states(ch, dch));
        std::vector<uint_t> M = dch.MMatrix();
        uint_t M_diag1 = dch.MDiag1();
        uint_t M_diag2 = dch.MDiag2();
        unsigned control_phase = ((M_diag1 >> 1) & one) + 2 * ((M_diag2 >> 1) & one);
        REQUIRE(control_phase == 3);
        REQUIRE(!!(M_diag1 & one));
    }
}

TEST_CASE("Random Deterministic Circuit")
{
    unsigned n_qubits = 10;
    CHState ch(n_qubits);
    DCHState dch(n_qubits);
    uint_t x_init = zer;
    std::string x_str = "";
    for(unsigned i=0; i<n_qubits; i++)
    {
        if(rand()%2)
        {
            x_init ^= (one << i);
            x_str += "1";
        }
        else
        {
            x_str += "0";
        }
    }
    CAPTURE(x_str);
    ch.CompBasisVector(x_init);
    dch.CompBasisVector(x_init);
    random_circuit_verify(ch, dch, deterministic_gates);
}

TEST_CASE("Test Hadamard Gates")
{
    unsigned n_qubits = 10;
    SECTION("Plus state")
    {
        CHState ch(n_qubits);
        DCHState dch(n_qubits);
        ch.H(0);
        dch.H(0);
        REQUIRE(check_states(ch, dch));
    }
    SECTION("Minus state")
    {
        CHState ch(n_qubits);
        DCHState dch(n_qubits);
        ch.X(0);
        dch.X(0);
        ch.H(0);
        dch.H(0);
        REQUIRE(check_states(ch, dch));
    }
}

TEST_CASE("Random H and Paulis")
{
    unsigned n_qubits = 10;
    CHState ch(n_qubits);
    DCHState dch(n_qubits);
    uint_t x_init = zer;
    std::string x_str = "";
    for(unsigned i=0; i<n_qubits; i++)
    {
        if(rand()%2)
        {
            x_init ^= (one << i);
            x_str += "1";
        }
        else
        {
            x_str += "0";
        }
    }
    // CAPTURE(x_str);
    std::vector<Gate> gs = {Gate::S, Gate::Sdag, Gate::H, Gate::X, Gate::Y, Gate::Z, Gate::CZ, Gate::H};
    random_circuit_verify(ch, dch, gs, 100);
}

scalar_t CHState::test_inner_product(uint_t &sample_1, uint_t &sample_2,
                                         std::vector<uint_t> &sample)
{
    scalar_t amp = this->InnerProduct(sample_1, sample_2, sample);
    return amp;
}

scalar_t DCHState::test_inner_product(uint_t &sample_1, uint_t &sample_2,
                                         std::vector<uint_t> &sample)
{
    scalar_t amp = this->InnerProduct(sample_1, sample_2, sample);
    return amp;
}

TEST_CASE("Inner Product Test")
{
    unsigned n_qubits = 10;
    CHState ch(10);
    DCHState dch(10);
    random_circuit_verify(ch, dch, all_gates, 30);
    REQUIRE(check_states(ch, dch));
    uint_t ds_1 = zer, ds_2 = zer;
    std::vector<uint_t> ds_ch(n_qubits, zer);
    std::vector<uint_t> ds_dch(n_qubits, zer);
    for(unsigned i=0; i<n_qubits; i++)
    {
        if(rand()%2)
        {
            ds_1 ^= (one << i);
            ds_ch[i] ^= (one << i);
        }
        if(rand()%2)
        {
            ds_2 ^= (one << i);
        }
        for(unsigned j=i+1; j<n_qubits; j++)
        {
            if(rand()%2)
            {
                ds_ch[i] ^= (one << j);
                ds_ch[j] ^= (one << i);
                ds_dch[i] ^= (one << j);
                ds_dch[j] ^= (one << i);
            }
        }
    }
    scalar_t ch_amp = ch.test_inner_product(ds_1, ds_2, ds_ch);
    ch_amp.Print();
    std::cout << std::endl;
    scalar_t dch_amp = dch.test_inner_product(ds_1, ds_2, ds_dch);
    dch_amp.Print();
    std::cout << std::endl;
    REQUIRE(std::abs(std::real(dch_amp.to_complex())-std::real(ch_amp.to_complex())) == Approx(0.));
    REQUIRE(std::abs(std::imag(dch_amp.to_complex())+std::imag(ch_amp.to_complex())) == Approx(0.));
}

int main( int argc, char* argv[] ) {
  time_t t;
  unsigned seed = (unsigned) time(&t);
  // unsigned seed = 1551981445;
  // unsigned seed = 1551995579;
  std::cout << "Initialized with seed: " << seed << std::endl;
  srand(seed);

  int result = Catch::Session().run( argc, argv );

  return result;
}

