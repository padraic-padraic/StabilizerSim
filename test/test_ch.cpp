#include "libstabilizer/core.hpp"
#include "libstabilizer/chstabilizer.hpp"
#include "libstabilizer/dchstabilizer.hpp"

#define CATCH_CONFIG_MAIN
#include "test/lib/catch.hpp"

#include <cstdlib>
#include <complex>
#include <vector>

namespace StabilizerSimulator
{

using complex_t = std::complex<double>;

TEST_CASE("Check computational amplitudes") 
{
    unsigned n_qubits = 5;
    uint_t dim = one << n_qubits;
    SECTION("Test initializer")
    {
        StabilizerState ch(n_qubits);
        DCHStabilizer dch(n_qubits);
        for(uint_t x=0; x<dim; x++)
        {
            scalar_t ch_amp = ch.Amplitude(x);
            scalar_t dch_amp = dch.Amplitude(x);
            if(ch_amp.eps != 0)
            {
                REQUIRE(ch_amp.p == dch_amp.p);
                REQUIRE(ch_amp.e == dch_amp.e);                
            }
            REQUIRE(ch_amp.eps == dch_amp.eps);
        }
    }
    SECTION("Test X")
    {
        StabilizerState ch(n_qubits);
        DCHStabilizer dch(n_qubits);
        for(unsigned i=0; i<n_qubits; i++)
        {
            if (rand() % 2)
            {
                ch.X(i);
                dch.X(i);
            }
        }
        for(uint_t x=0; x<dim; x++)
        {
            scalar_t ch_amp = ch.Amplitude(x);
            scalar_t dch_amp = dch.Amplitude(x);
            if(ch_amp.eps != 0)
            {
                REQUIRE(ch_amp.p == dch_amp.p);
                REQUIRE(ch_amp.e == dch_amp.e);                
            }
            REQUIRE(ch_amp.eps == dch_amp.eps);
        }
    }
}

}