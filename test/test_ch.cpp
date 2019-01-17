#include "libstabilizer/core.hpp"
#include "libstabilizer/chstabilizer.hpp"
#include "libstabilizer/dchstabilizer.hpp"

#define CATCH_CONFIG_RUNNER
#include "test/lib/catch.hpp"

#include <algorithm>
#include <cstdlib>
#include <complex>
#include <iostream>
#include <time.h>
#include <vector>

using namespace StabilizerSimulator;

using complex_t = std::complex<double>;

const double precision=1e-8;

bool complex_close(complex_t &a, complex_t &b)
{
    return( (std::abs(a.real()-b.real())<1e-8) & (std::abs(a.imag()-b.imag())<1e-8) );
}


bool check_states(StabilizerState &ch, DCHStabilizer &dch)
{
    unsigned n = ch.NQubits();
    uint_t dim = one << n;
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
            return false;
        }
    }
    return true;
}

TEST_CASE("Test initialisation and Amplitude")
{
    unsigned n_qubits = 5;
    StabilizerState ch(n_qubits);
    DCHStabilizer dch(n_qubits);
    REQUIRE(check_states(ch, dch));
}

TEST_CASE("Test Computational state initialisation")
{
    unsigned n_qubits = 10;
    StabilizerState ch(n_qubits);
    DCHStabilizer dch(n_qubits);
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

TEST_CASE("Test Z phases")
{
    unsigned n_qubits = 10;
    StabilizerState ch(n_qubits);
    DCHStabilizer dch(n_qubits);
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

TEST_CASE("Test Pauli Y")
{
    unsigned n_qubits = 10;
    StabilizerState ch(n_qubits);
    DCHStabilizer dch(n_qubits);
    for(unsigned i=0; i<20; i++)
    {
        unsigned target = (rand() % 10);
        INFO("Target is " << target);
        ch.Y(target);
        dch.Y(target);
        REQUIRE(check_states(ch, dch));
    }
}

TEST_CASE("Test S Gates")
{
    unsigned n_qubits = 10;
    StabilizerState ch(n_qubits);
    DCHStabilizer dch(n_qubits);
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
    StabilizerState ch(n_qubits);
    DCHStabilizer dch(n_qubits);
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

TEST_CASE("Test CX gate")
{
    unsigned n_qubits = 10;
    SECTION("Check computational strings")
    {
        StabilizerState ch(n_qubits);
        DCHStabilizer dch(n_qubits);
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
        REQUIRE(check_states(ch, dch));
        // Check trivial CX
        ch.CX(8, control);
        dch.CX(8, control);
        REQUIRE(check_states(ch, dch));
    }
    SECTION("Check phases")
    {
        StabilizerState ch(n_qubits);
        DCHStabilizer dch(n_qubits);
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
        REQUIRE(!!(M[0] & (one << 1)));
        REQUIRE(!!(M_diag1 & one));
    }
}

TEST_CASE("Test Hadamard Gates")
{
    unsigned n_qubits = 10;
    StabilizerState ch(n_qubits);
    DCHStabilizer dch(n_qubits);
    ch.H(0);
    dch.H(0);
    REQUIRE(check_states(ch, dch));
}

int main( int argc, char* argv[] ) {
  time_t t;
  unsigned seed = (unsigned) time(&t);
  std::cout << "Initialized with seed: " << seed << std::endl;
  srand(seed);

  int result = Catch::Session().run( argc, argv );

  return result;
}

