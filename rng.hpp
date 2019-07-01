#ifndef RNG_HPP
#define RNG_HPP
#include <stdint.h>

//===============================//
//          RNG Code             //
//===============================//

/*  Written in 2018 by David Blackman and Sebastiano Vigna (vigna@acm.org)

To the extent possible under law, the author has dedicated all copyright
and related and neighboring rights to this software to the public domain
worldwide. This software is distributed without any warranty.

See <http://creativecommons.org/publicdomain/zero/1.0/>. */

/* This is a fixed-increment version of Java 8's SplittableRandom generator
   See http://dx.doi.org/10.1145/2714064.2660195 and 
   http://docs.oracle.com/javase/8/docs/api/java/util/SplittableRandom.html

   It is a very fast generator passing BigCrush, and it can be useful if
   for some reason you absolutely want 64 bits of state; otherwise, we
   rather suggest to use a xoroshiro128+ (for moderately parallel
   computations) or xorshift1024* (for massively parallel computations)
   generator. */
namespace StabilizerSimulator
{

thread_local static uint64_t sm_x; /* The state can be seeded with any value. */

uint64_t sm_next() {
  uint64_t z = (sm_x += 0x9e3779b97f4a7c15);
  z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9;
  z = (z ^ (z >> 27)) * 0x94d049bb133111eb;
  return z ^ (z >> 31);
}

/* This is xoshiro256** 1.0, our all-purpose, rock-solid generator. It has
   excellent (sub-ns) speed, a state (256 bits) that is large enough for
   any parallel application, and it passes all tests we are aware of.

   For generating just floating-point numbers, xoshiro256+ is even faster.

   The state must be seeded so that it is not everywhere zero. If you have
   a 64-bit seed, we suggest to seed a splitmix64 generator and use its
   output to fill s. */
static inline uint64_t xoro_rotl(const uint64_t x, int k) {
    return (x << k) | (x >> (64 - k));
}

thread_local static uint64_t xoro_s[4];

uint64_t xoro_next(void) {
    const uint64_t result_starstar = xoro_rotl(xoro_s[1] * 5, 7) * 9;

    const uint64_t t = xoro_s[1] << 17;

    xoro_s[2] ^= xoro_s[0];
    xoro_s[3] ^= xoro_s[1];
    xoro_s[1] ^= xoro_s[2];
    xoro_s[0] ^= xoro_s[3];

    xoro_s[2] ^= t;

    xoro_s[3] = xoro_rotl(xoro_s[3], 45);

    return result_starstar;
}

//===============================//
//        Wrapper funcs          //
//===============================//

void init_rng(unsigned seed, unsigned offset=0)
{
  sm_x = (seed + offset);
  xoro_s[0] = sm_next();
  xoro_s[1] = sm_next();
  xoro_s[2] = sm_next();
  xoro_s[3] = sm_next();
}

double random_double()
{
  uint64_t x = xoro_next();
  return (x >> 11) * (1. / (UINT64_C(1) << 53));
}

bool random_bit()
{
  return (xoro_next() & (1ULL));
}

uint64_t random_uint()
{
  return xoro_next();
}
}
#endif