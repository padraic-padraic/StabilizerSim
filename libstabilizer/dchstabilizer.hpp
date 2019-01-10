#ifndef DCH_STABILIZER_HPP
#define DCH_STABILIZER_HPP

#include "Core.hpp"

#include <algorithm>
#include <vector>

const char* dch_fnames[] = {"s", "v", "M", "M_diag", "G"};

namespace StabilizerSimlator
{

class DCHStabilizer
{
public:
    // Construct DCH Stabilizer corresponding to |0>^{\otimes n}
    DCHStabilizer(unsigned n_qubits);
    //Copy constructor
    DCHStabilizer(const DCHStabilizer& rhs);

     uint_t NQubits() const
      {
        return n;
      }
      scalar_t Omega() const
      {
        return omega;
      }

      uint_t M_diag1() const
      {
        return M_diag1;
      }

      uint_t M_diaga2() const
      {
        return M_diag2;
      }

      std::vector<uint_t> GMatrix() const
      {
        return G;
      }

      std::vector<uint_t> MMatrix() const
      {
        return M;
      }

    void CX(unsigned control, unsigned target);
    void CZ(unsigned control, unsigned target);
    void S(unsigned target);
    void Sdag(unsigned target);
    void Z(unsigned target);
    void H(unsigned target);

    // state preps
    void CompBasisVector(uint_t x); // prepares |x> 
    void HadamardBasisVector(uint_t x); // prepares H(x)|00...0>

    //Measurements
    scalar_t Amplitude(uint_t x);

private:
    unsigned n; //N Qubits
    uint_t s; //Computational basis string |s>
    uint_t v; //Hadamard string: U_{H}|s> = H(v)|s>
    std::vector<uint_t> G; // CNOT matrix G: U_{c} |s> = |sG>
    std::vector<uint_t> G_inverse;
    //Store the phase matrix M: U_{D}|x> = i^{xDx^{T}}|x>
    uint_t M_diag1;// Diagonal elements of the phase matrix M, mod 4.
    uint_t M_diag2;
    std::vector<uint_t> M;// CZ part of the phase matrix M
    scalar_t omega; //Internal phase: |phi> = omega*U_{d}*U_{c}*U_{h}|s>

    std::vector<uint_t> GT;
    std::vector<uint_t> G_inverseT;
    bool isReadyGT;
    bool isReadyG_inverseT;
    void TransposeG();
    void TransposeG_inverse();
    //Commute a Pauli through the D, C and H layers
    void CommutePauli(pauli_t& P);
    pauli_t GetAmplitudeX(uint_t x);
    //Update used for Hadamards and Measurement
    void UpdateSVector(uint_t t, uint_t u, unsigned b);
};

DCHStabilizer::DCHStabilizer(unsigned n_qubits) :
n(n_qubits),
s(zer),
v(zer),
G(n_qubits,zer),
M_diag1(zer),
M_diag2(zer),
M(n, zer),
GT(n_qubits,zer),
G_inverseT(n_qubits, zer),
isReadyGT(false),
isReadyG_inverseT(false)
{
  // Initialize G to the identity
  for (unsigned i=0; i<n_qubits; i++)
  {
    G[i] |= (one << i);
  }
};

DCHStabilizer::DCHStabilizer(const DCHStabilizer& rhs) :
n(rhs.n),
s(rhs.s),
v(rhs.v),
G(rhs.G),
G_inverse(rhs.G_inverse);
M_diag1(rhs.M_diag1),
M_diag2(rhs.M_diag2),
M(rhs.M),
GT(rhs.GT),
G_inverseT(rhs.G_inverseT),
isReadyGT(rhs.isReadyGT),
isReadyG_inverseT(rhs.isReadyG_inverseT)
{
};

void DCHStabilizer::CompBasisVector(uint_t x)
{
  s = x;
  M_diag1 = zer;
  M_diag2 = zer;
  std::fill(M.begin(), M.end(), zer);
  for (unsigned i=0; i<n; i++)
  {
    G[i] = (one << i);
    G_inverse[i] = (one << i);
  }
  isReadyGT = false;
  isReadyG_inverseT = false;
}

void DCHStabilizer::HadamardBasisVector(uint_t x)
{
  s = zer;
  v = x;
  M_diag1 = zer;
  M_diag2 = zer;
  std::fill(M.begin(), M.end(), zer);
  for (unsigned i=0; i<n; i++)
  {
    G[i] = (one << i);
    G_inverse[i] = (one << i);
  }
  isReadyGT = false;
  isReadyG_inverseT = false;
}

void DCHStabilizer::S(unsigned target)
{
  // Add one to M_{ii}
  if ((M_diag1 >> target) & one)
  {
    M_diag2 ^= (one << target);
  }
  M_diag1 ^= (one << target);
}

void DCHStabilizer::Z(unsigned target)
{
  //Add two to M_{ii}
  M_diag2 ^= (one << target);
}

void DCHStabilizer::Sdag(unsigned target)
{
  //Add three to M_{ii}
  // 0 + 3 = 3; => 00 + 11 = 11;
  // 1 + 3 = 4 = 0; => 01 + 11 = 00;
  // 2 + 3 = 5 = 1; => 10 + 11 = 01;
  // 3 + 3 = 6 = 2; => 11 + 11 = 10;
  // Simply flip the ones bit, conditionally flip twos bit
  if (!((M_diag1 >> target) & one))
  {
    M_diag2 ^= (one << target);
  }
  M_diag1 ^= (one << target);
}

void DCHStabilizer::CZ(unsigned control, unsigned target)
{
  M[target] ^= (one << control);
  M[control] ^= (one << target);
}

void DCHStabilizer::CX(unsigned control, unsigned target)
{

  isReadyGT = false;
  uint_t target_col = M[target];
  uint_t shift = (one << control);
  //Update the con
  M[control] ^= target_col;
  for (unsigned i=0; i<n; i++)
  {
    // Update the control row of M
    M[i] ^= (((target_col >> i) & one) * shift);
    //Update the control row of G^-1
    G_inverse[i] ^= ( ((G_inverse[i] >> target) & one)  * shift);
  }
  // M_{c,c} += 2* M_{c,t}
  if ((target_col >> control) & one)
  {
    M_diag2 ^= shift;
  }
  // M_{c,c} += M_{t,t};
  uint_t diagonal_target_bit = ((M_diag1 >> target) & one);
  if ((M_diag1 >> control) & diagonal_target_bit)
  {
    M_diag2 ^= shift;
  }
  M_diag1 ^= diagonal_target_bit * shift;
  M[control] ^= diagonal_target_bit * shift;
  M_diag2 ^= ((M_diag2 >> target) & one) * shift;
  // Update the target column of G
  G[target] ^= G[control];
}

void DCHStabilizer::TransposeG()
{
  for (unsigned i=0; i<n; i++)
  {
    uint_t shift = (one << i);
    for (unsigned j=0; j<n; j++)
    {
      if((G[i] >>j) & one)
      {
        GT[j] |= shift;
      } else {
        GT[j] &= ~(shift);
      }
    }
  }
  isReadyGT = true;
}

void DCHStabilizer::TransposeG_inverse()
{
  for (unsigned i=0; i<n; i++)
  {
    uint_t shift = (one << i);
    for (unsigned j=0; j<n; j++)
    {
      if((G_inverse[i] >>j) & one)
      {
        G_inverseT[j] |= shift;
      } else {
        G_inverseT[j] &= ~(shift);
      }
    }
  }
  isReadyG_inverseT = true;
}

void DCHStabilizer::CommutePauli(pauli_t& p)
{
  uint_t x_temp = zer;
  uint_t z_temp = zer;
  char phase = hamming_weight(p.X&M_diag1) + 2*hamming_weight(p.X&M_diag2);
  for (unsigned i=0; i<n; i++)
  {
    p.Z ^= (one << i) * hamming_parity(p.X&M[i]);
  }
  //Shift phase to mod 8
  phase = (phase % 4) *2;
  p.e -= phase;
  p.e = p.e % 8;
  if (p.e < 0)
  {
    p.e += 8;
  }
  //Commute through C:
  for (unsigned i=0; i<n; i++)
  {
    uint_t shift = (one << i);
    if (!!(p.Z&shift))
    {
      ztemp ^= G[i];
    }
    xtemp ^= shift * (hamming_weight(p.X&G_inverse[i]) & one);
  }
  //Commute through H
  p.X = (ztemp & v) ^ (xtemp & ~v);
  p.Z = (xtemp & v) ^ (ztemp & ~v);
}

pauli_t DCHStabilizer::GetAmplitudeX(uint_t x) {
  pauli_t p;
  p.X = x;
  char phase = hamming_weight(p.x&M_diag1) + 2*hamming_weight(p.x&M_diag2);
  for (unsigned i=0; i<n; i++)
  {
    p.Z ^= (one << i) * hamming_parity(p.X&M[i]);
  }
  //Shift phase to mod 8
  phase = (phase % 4) *2;
  p.e -= phase;
  p.e = p.e % 8;
  if (p.e < 0)
  {
    p.e += 8;
  }
  //Commute through C:
  uint_t x_temp = zer;
  uint_t z_temp = zer;
  for (unsigned i=0; i<n; i++)
  {
    uint_t shift = (one << i);
    if (!!(p.Z&shift))
    {
      ztemp ^= G[i];
    }
    xtemp ^= shift * (hamming_weight(p.X&G_inverse[i]) & one);
  }
  p.X = x_temp;
  p.Z = z_temp;
  return p;
}

scalar_t DCHStabilizer::Amplitude(uint_t x)
{
  scalar_t amp;
  if(omega.eps == 0)
  {
    amp.eps = 0;
    return amp;
  }
  amp.p = -1 * ((int) hamming_weight(v));
  pauli_t p = GetAmplitudeX(x);
  amp.e = 2*p.e;
  if (!!((p.X ^ s) & ~v))
  {
    amp.eps = 0;
    return amp;
  }
  if (hamming_parity(p.X&s&v))
  {
    amp.e += 4;
    amp.e = (amp.e%8);
  }
  amp.conjugate();
  amp *= omega;
  return amp;
}

void DCHStabilizer::UpdateSVector(uint_t t, uint_t u, unsigned b)
{
  b %=4;
  if (t==u) {
    switch(b)
    {
      case 0:
        omega.p +=1;
        s=t;
        return;
      case 1:
        s=t;
        omega.e = (omega.e+1)%8;
        return;
      case 2:
        s=t;
        omega.eps = 0;
        return;
      case 3:
        s=t;
        omega.e = (omega.e +  7) % 8;
        return;
      default:
        assert(0);
    }
  }
  uint_t nu0 = (t^u) & (~v);
  uint_t nu1 = (t^u) & v;
  //Pick 
  unsigned q = 0;
  uint_t q_shift = zer;
  bool nu0Empty = true;
  if (nu0)
  {
    nu0Empty = false;
    //nu0 is non-empty, we need to build nu1
    //and then commute the pseudoCZ gate through the
    //current Uc

    //First up, pick q
    for (unsigned i=0; i<n; i++)
    {
      if (!!(nu0 & (one << q)))
      {
        q=i;
        q_shift = (one << q);
        break;
      }
    }
    //
    assert((nu0 & q_shift)>0);
    nu0 ^= q_shift;
  } else {
    //nu0 is empty
    // Pick q
    for (unsigned i=0; i<n; i++)
    {
      if (!!(nu1 & (one << q)))
      {
        q=i;
        q_shift = (one << q);
        break;
      }
    }
    assert((nu1 & q_shift)>0);
    nu1 ^= q_shift;
  }

  //
  if (t & q_shift)
    {
        s=u;
        omega.e=(omega.e + 2*b) % 8;
        b=(4-b) % 4;
        assert(!(u & q_shift));
    }
    else
    {
        s=t;
    }

    // change the order of H and S gates
    // H^{a} S^{b} |+> = eta^{e1} S^{e2} H^{e3} |e4>
    // here eta=exp( i (pi/4) )
    // a=0,1
    // b=0,1,2,3
    //
    // H^0 S^0 |+> = eta^0 S^0 H^1 |0>
    // H^0 S^1 |+> = eta^0 S^1 H^1 |0>
    // H^0 S^2 |+> = eta^0 S^0 H^1 |1>
    // H^0 S^3 |+> = eta^0 S^1 H^1 |1>
    //
    // H^1 S^0 |+> = eta^0 S^0 H^0 |0>
    // H^1 S^1 |+> = eta^1 S^1 H^1 |1>
    // H^1 S^2 |+> = eta^0 S^0 H^0 |1>
    // H^1 S^3 |+> = eta^{7} S^1 H^1 |0>
    //
    // "analytic" formula:
    // e1 = a * (b mod 2) * ( 3*b -2 )
    // e2 = b mod 2
    // e3 = not(a) + a * (b mod 2)
    // e4 = not(a)*(b>=2) + a*( (b==1) || (b==2) )

    bool a=((v & q_shift)>0);
    unsigned e1=a*(b % 2)*( 3*b -2);
    unsigned e2 = b % 2;
    bool e3 = ( (!a) != (a && ((b % 2)>0) ) );
    bool e4 = ( ( (!a) && (b>=2) ) != (a && ((b==1) || (b==2)))  );
    //Split up based on if nu0 is empty
    if (nu0Empty)
    {
      if (nu1)//if T2 != Identity
      {
        //Right multiply UC by T2
        for (unsigned i=0; i<n; i++)
        {
          if((nu1 >> i) & one)
          {
            G[q] ^= G[i];
          }
        }
        if(e2)
        {
          uint_t rowQ = G[q];
          // We need to commute S through and then add it to UD
          for (unsigned i=0; i<n; i++)
          {
            uint_t shift = (one << i);
            if (!!(rowQ & shift))
            {
              if (!!(M_diag1 & shift))
              {
                M_diag2 ^= shift;
              }
              M_diag1 ^= shift;
              for (unsigned j=i+1; j<n; j++)
              {
                if ((rowQ >> j) & one)
                {
                  M[i] ^= (one << j);
                  M[j] ^= shift;
                }
              }
            }
          }
        }
      } else {
        //T2 is identity, but we still need to add S to UD
        if (e2)
        {
          M_diag2 ^= (M_diag1 & (one << q));
          M_diag1 ^= (one << q);
        }
      }
    } else {
      //nu0 Not Empty
      if (nu0)
      {
        //Commute WD through UC and add to UD
        //Setup the Z vector, pull out the y vector
        uint_t y = G[q];
        uint_t z = zer;
        for (unsigned j=0; j<n; j++)
        {
          if((nu1>>j) & one)
          {
            z ^= G[j];
          }
        }
        // Diagonal elemenets M_{ii} += y_{i} = 2*z_{i}y_{i}
        // 2* z_{i} = 0 if Z=0,2 or  2 if Z=1,3
        if (e2)
        {
          M_diag1 ^= y;          
        }
        M_diag2 ^= (z & y);
        // Offdiagonal elements
        for (unsigned i=0; i<n; i++)
        {
          bool y_i = (y>>i)&one;
          bool z_i = (z>>i)&one;
          uint_t shift = one << i;
          for (unsigned j=i+1; j<n; j++)
          {
            bool y_j = (y>>j)&one;
            bool z_j = (z>>j)&one;
            bool val = (e2 & (y_i & y_j));
            val ^= ((y_i&z_j) ^ (z_i &y_j));
            M[j] ^= val*shift;
            M[i] ^= val*(one << j);
          }
        }
      }
    }
    // set q-th bit of s to e4
    s&=~q_shift;
    s^=e4*q_shift;

    // set q-th bit of v to e3
    v&=~q_shift;
    v^=e3*q_shift;

    // update the scalar factor omega
    omega.e=(omega.e  + e1) % 8;
}

}
#endif