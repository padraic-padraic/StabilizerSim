#ifndef DCH_STABILIZER_HPP
#define DCH_STABILIZER_HPP

#include "core.hpp"

#include <algorithm>
#include <vector>

const char* dch_fnames[] = {"s", "v", "M", "M_diag", "G"};

namespace StabilizerSimulator
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

      uint_t MDiag1() const
      {
        return M_diag1;
      }

      uint_t MDiag2() const
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
    void X(unsigned target);
    void Y(unsigned target);
    void H(unsigned target);

    // state preps
    void CompBasisVector(uint_t x); // prepares |x> 
    void HadamardBasisVector(uint_t x); // prepares H(x)|00...0>

    //Measurements
    scalar_t Amplitude(uint_t x);
    void MeasurePauli(pauli_t P); // applies a gate (I+P)/2 
                                // where P is an arbitrary Pauli operator
    void MeasurePauliProjector(std::vector<pauli_t>& generators);

    #ifdef CATCH_VERSION_MAJOR //Helper 
    void test_commute_pauli();
    #endif

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
    bool isReadyGT;
    void TransposeG();
    //Commute a Pauli through the D, C and H layers
    void CommutePauliDC(pauli_t &P);
    void CommutePauli(pauli_t &P);
    //Update used for Hadamards and Measurement
    void UpdateSVector(uint_t t, uint_t u, unsigned b);
};

//-------------------------------//
// Implementation                //
//-------------------------------//

DCHStabilizer::DCHStabilizer(unsigned n_qubits) :
n(n_qubits),
s(zer),
v(zer),
G(n_qubits,zer),
G_inverse(n_qubits, zer),
M_diag1(zer),
M_diag2(zer),
M(n, zer),
GT(n_qubits,zer),
isReadyGT(false)
{
  // Initialize G to the identity
  for (unsigned i=0; i<n_qubits; i++)
  {
    G[i] |= (one << i);
    G_inverse[i] |= (one << i);
  }
};

DCHStabilizer::DCHStabilizer(const DCHStabilizer& rhs) :
n(rhs.n),
s(rhs.s),
v(rhs.v),
G(rhs.G),
G_inverse(rhs.G_inverse),
M_diag1(rhs.M_diag1),
M_diag2(rhs.M_diag2),
M(rhs.M),
GT(rhs.GT),
isReadyGT(rhs.isReadyGT)
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

void DCHStabilizer::Z(unsigned target)
{
  //Add two to M_{ii}
  M_diag2 ^= (one << target);
}

void DCHStabilizer::X(unsigned target)
{
  pauli_t p;
  p.X ^= (one << target);
  CommutePauli(p);
  omega.e += (hamming_parity(p.Z&s))*4;
  omega.e += 2*(p.e);
  omega.e = omega.e % 8;
  s ^= p.X;
}

void DCHStabilizer::Y(unsigned target)
{
  Z(target);
  X(target);
  omega.e = (omega.e + 2) % 8;
}

void DCHStabilizer::H(unsigned target)
{
  // Represent H_{a} = \frac{1}{2}\left(X_{a} + Z_{a}\right)
  std::cout << "Hadamard on qubit " << target << std::endl;
  pauli_t p, q;
  p.X ^= (one << target);
  q.Z ^= (one << target);
  CommutePauli(p);
  CommutePauli(q);
  uint_t t = s^p.X;
  unsigned t_phase = hamming_parity(s&p.Z) * 2U;
  t_phase += p.e;
  t_phase = t_phase % 4;
  uint_t u = s^q.X;
  unsigned u_phase = hamming_parity(s&q.Z) * 2U;
  if(q.e != 0)
  {
    throw std::logic_error("This pauli cannot have non-zero phase.");
  }
  unsigned b;
  if(u_phase == 0)
  {
    b = (4-t_phase);
  }
  else
  {
    switch(t_phase)
    {
      case 0:
        b = 2;
        break;
      case 1:
        b = 1;
        break;
      case 2:
        b=0;
        break;
      case 3:
        b=3;
        break;
      default:
        throw std::logic_error("Wat");
    }
  }
  std::cout << "Global phase is: " << t_phase << std::endl;
  omega.e = (omega.e + 2*t_phase) %8;
  b = b%4;
  if(t==u)
  {
    s = t;
    if(!((b==1) || (b==3))) // otherwise the state is not normalized
    {
      throw std::logic_error("State is not properly normalised, b should be 1 or 3.\n");
    }
    if (b==1)
    {
        omega.e=(omega.e + 1) % 8;
    }
    else
    {
        omega.e=(omega.e + 7) % 8;
    }
  }
  else
  {
    UpdateSVector(t, u, b);
  }
}

void DCHStabilizer::CZ(unsigned control, unsigned target)
{
  if(control == target)
  {
    throw std::logic_error("Controlled operation cannot target the same qubit.");
  }
  M[target] ^= (one << control);
  M[control] ^= (one << target);
}

void DCHStabilizer::CX(unsigned control, unsigned target)
{
  if(control == target)
  {
    throw std::logic_error("Controlled operation cannot target the same qubit.");
  }
  isReadyGT = false;
  uint_t target_col = M[target];
  uint_t shift = (one << control);
  //Update the control col of M
  M[control] ^= target_col;
  bool target_bit = !!(M_diag1 & (one << target));
  for (unsigned i=0; i<n; i++)
  {
    // Update the control row of M
    // We use the target col as M is symmetric
    M[i] ^= (((target_col >> i) & one) * shift);
    //Update the control row of G^-1
    G_inverse[i] ^= ( ((G_inverse[i] >> target) & one)  * shift);
  }
  //Update control row/column with M_{t,t}
  M[control] ^= target_bit * (one << target);
  M[target] ^= target_bit * shift;
  // M_{c,c} += 2* M_{c,t}
  if ((target_col >> control) & one)
  {
    M_diag2 ^= shift;
  }
  // M_{c,c} += M_{t,t};
  if(!!(M_diag1 & shift) & target_bit)
  {
    M_diag2 ^= shift;
  }
  M_diag1 ^= target_bit*shift;
  M_diag2 ^= ((M_diag2 >> target) & one)*shift;
  // Update the target column of G
  G[target] ^= G[control];
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
  pauli_t p;
  p.X = x;
  CommutePauliDC(p);
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

void DCHStabilizer::MeasurePauli(pauli_t P)
{
  CommutePauli(P);
  unsigned b = P.e;
  uint_t u = s^P.X;
  b += hamming_parity(s&P.Z) * 2;
  b = b%4;
  UpdateSVector(s, u, b);
}

void DCHStabilizer::MeasurePauliProjector(std::vector<pauli_t>& generators)
{
  for (uint_t i=0; i<generators.size(); i++)
  {
      this->MeasurePauli(generators[i]);
      if (omega.eps == 0)
      {
          break;
      }
  }
}

//----------------------------------//
// Implementation - private methods //
//----------------------------------//

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

void DCHStabilizer::CommutePauliDC(pauli_t &p)
{
  unsigned phase = hamming_weight(p.X&M_diag1) + 2*hamming_weight(p.X&M_diag2);
  uint_t offdiag_terms = zer;
  for (unsigned i=0; i<n; i++)
  {
    p.Z ^= (one << i) * hamming_parity(p.X&M[i]);
    if((p.X>>i)& one)
    {
      offdiag_terms += hamming_weight(p.X&M[i]);
    }
  }
  //Compensate for the zero diagonal in M
  p.Z ^= (p.X & M_diag1);
  phase += (!!(offdiag_terms & (one << 1))) * 2U;
  phase = (phase % 4);
  phase = (4-phase) %4;
  p.e = phase;
  //Commute through C:
  uint_t x_temp = zer;
  uint_t z_temp = zer;
  for (unsigned i=0; i<n; i++)
  {
    uint_t shift = (one << i);
    if (!!(p.Z&shift))
    {
      z_temp ^= G[i];
    }
    x_temp ^= shift * hamming_parity(p.X & G_inverse[i]);
  }
  p.X = x_temp;
  p.Z = z_temp;
}

void DCHStabilizer::CommutePauli(pauli_t& p)
{
  CommutePauliDC(p);
  uint_t x_temp = p.X;
  uint_t z_temp = p.Z;
  //Commute through H
  p.X = (z_temp & v) ^ (x_temp & ~v);
  p.Z = (x_temp & v) ^ (z_temp & ~v);
}

void DCHStabilizer::UpdateSVector(uint_t t, uint_t u, unsigned b)
{
  b %= 4;
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
        throw std::logic_error("Invalid phase factor found b:" + std::to_string(b) + ".\n");
    }
  }
  std::cout << "T string is :";
  Print(t, n);
  std::cout << std::endl;
  std::cout << "U string is :";
  Print(u, n);
  std::cout << std::endl;
  uint_t nu0 = (t^u) & (~v);
  uint_t nu1 = (t^u) & v;
  unsigned q = 0;
  uint_t q_shift = zer;
  bool nu0Empty = true;
  std::cout << "nu0: ";
  Print(nu0, n);
  std::cout << std::endl;
  std::cout << "nu1: ";
  Print(nu1, n);
  std::cout << std::endl;
  if (nu0)
  {
    nu0Empty = false;
    //nu0 is non-empty, we need to build nu1
    //and then commute the pseudoCZ gate through the
    //current Uc

    //First up, pick q
    for (unsigned i=0; i<n; i++)
    {
      if (!!(nu0 & (one << i)))
      {
        q=i;
        q_shift = (one << i);
        break;
      }
    }
    //
    nu0 ^= q_shift;
  }
  else
  {
    //nu0 is empty
    // Pick q
    for (unsigned i=0; i<n; i++)
    {
      if (!!(nu1 & (one << i)))
      {
        q=i;
        q_shift = (one << i);
        break;
      }
    }
    nu1 ^= q_shift;
  }
  //Pick the string with the q-th bit equal to zero. If necessary, swap the strings.
  if (!!(t & q_shift))
  {
    std::cout << "Swapping strings based on the qth bit" << std::endl;
    s=u;
    b=(4-b) % 4;
    if(!!(u & q_shift))
    {
      throw std::logic_error("t and u strings do not differ.");
    }
  }
  else
  {
      s=t;
  }
  //We know have the qth bit looks like 
  // |0>+i^{b}|1> = S^{b}|+>

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
  std::cout << "H/S params: e1 "<< e1 << " e2:" << e2 << " e3:" << e3 << " e4:" << e4 << std::endl;
  //Split up based on if nu0 is empty
  if (nu0Empty)
  {
    std::cout << "nu0 is empty" << std::endl;
    if (nu1)//if T2 != Identity
    {
      //Right multiply UC by T2
      for(unsigned i=0; i<n; i++)
      {
        if((nu1 >> i) & one)
        {
          //xor col i of G^{-1} into col q
          G_inverse[q] ^= G_inverse[i];
          uint_t shift = (one << i);
          for(unsigned j=0; j<n; j++)
          {
            //xor row q of G into row i
            G[j] ^= ((!!(G[j] & q_shift)) * shift);
          }
        }
      }
    }
    if(e2) // Commute the s gate through Uc 
    {
      uint_t y = G_inverse[q];
      M_diag2 ^= (y & M_diag1);
      M_diag1 ^= y;
      for(unsigned i=0; i<n; i++)
      {
        M[i] ^= ((y>>i)&one)*y;
        M[i] &= ~(one << i);
      }
    }
  }
  else
  {
    std::cout << "nu0 not empty" << std::endl;
    //Commute WD through UC and add to UD
    //Setup the Z vector, pull out the y vector
    uint_t y = G_inverse[q];
    uint_t z = zer;
    for (unsigned j=0; j<n; j++)
    {
      if((nu1>>j) & one)
      {
        z ^= G_inverse[j];
      }
    }
    // Diagonal elemenets M_{ii} += y_{i} + 2*z_{i}y_{i}
    // 2* z_{i} = 0 if Z=0,2 or  2 if Z=1,3
    if (e2)
    {
      M_diag2 ^= (y&M_diag1);
      M_diag1 ^= y;
    }
    M_diag2 ^= (z & y);
    // Offdiagonal elements
    for (unsigned i=0; i<n; i++)
    {
      if((y>>i)&one)
      {
        M[i] ^= z;
        if(e2)
        {
          M[i] ^= y;
        }
      }
      if((z>>i)&one)
      {
        M[i] ^= y;
      }
      M[i] &= ~(one << i);
    }
    //Right multiply G by T0
    if(nu0)
    {
      for(unsigned i=0; i<n; i++)
      {
        if((nu0 >> i) & one)
        {
          //xor col i of G^{-1} into col q
          G_inverse[q] ^= G_inverse[i];
          uint_t shift = (one << i);
          for(unsigned j=0; j<n; j++)
          {
            //xor row q of G into row i
            G[j] ^= ((!!(G[j] & q_shift)) * shift);
          }
        }
      }
    }
  }
    
  // set q-th bit of s to e4
  s&=~q_shift;
  s^= (e4*q_shift);
  // set q-th bit of v to e3
  v&=~q_shift;
  v^=e3*q_shift;
  std::cout << "V after:";
  Print(v, n);
  std::cout << std::endl;

  // update the scalar factor omega
  omega.e=(omega.e  + e1) % 8;
}

}
#endif