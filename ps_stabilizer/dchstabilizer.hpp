#ifndef DCH_STABILIZER_HPP
#define DCH_STABILIZER_HPP

#include "core.hpp"

#include <vector>

const char* dch_fnames[] = {"s", "v", "M", "M_diag", "G"};

namespace StabilizerSimulator
{

class DCHState
{
public:
    // Construct DCH Stabilizer corresponding to |0>^{\otimes n}
  DCHState(unsigned n_qubits);
  //Copy constructor
  DCHState(const DCHState& rhs);

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

  scalar_t InnerProduct(DCHState& other);

  template<class T> friend double NormEstimate(std::vector<T>& states,
          const std::vector< std::complex<double> >& phases, 
          const std::vector<uint_t>& Samples_d1,
          const std::vector<uint_t> &Samples_d2, 
          const std::vector< std::vector<uint_t> >& Samples);
  template<class T> friend double NormEstimate(std::vector<T>& states,
          const std::vector< std::complex<double> >& phases, 
          const std::vector<uint_t>& Samples_d1,
          const std::vector<uint_t> &Samples_d2, 
          const std::vector< std::vector<uint_t> >& Samples,
          int n_threads);

  #ifdef CATCH_VERSION_MAJOR //Helper 
  void test_commute_pauli();
  scalar_t test_inner_product(uint_t &sample_1, uint_t &sample_2,
          std::vector<uint_t> &sample);
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

  scalar_t InnerProduct(uint_t &A_diag1, uint_t &A_diag2, std::vector<uint_t> &A);
};

//-------------------------------//
// Implementation                //
//-------------------------------//

DCHState::DCHState(unsigned n_qubits) :
n(n_qubits),
s(zer),
v(zer),
G(n_qubits,zer),
G_inverse(n_qubits, zer),
M_diag1(zer),
M_diag2(zer),
M(n, zer),
GT(n_qubits, zer),
isReadyGT(true)
{
  // Initialize G to the identity
  for (unsigned i=0; i<n_qubits; i++)
  {
    G[i] |= (one << i);
    GT[i] |= (one << i);
    G_inverse[i] |= (one << i);
  }
};

DCHState::DCHState(const DCHState& rhs) :
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

void DCHState::CompBasisVector(uint_t x)
{
  s = x;
  M_diag1 = zer;
  M_diag2 = zer;
  for (unsigned i=0; i<n; i++)
  {
    M[i] = zer;
    G[i] = (one << i);
    GT[i] = (one << i);
    G_inverse[i] = (one << i);
  }
  isReadyGT = true;
}

void DCHState::HadamardBasisVector(uint_t x)
{
  s = zer;
  v = x;
  M_diag1 = zer;
  M_diag2 = zer;
  for (unsigned i=0; i<n; i++)
  {
    M[i] = zer;
    G[i] = (one << i);
    G_inverse[i] = (one << i);
  }
  isReadyGT = false;
}

void DCHState::S(unsigned target)
{
  // Add one to M_{ii}
  uint_t shift = (one << target);
  M_diag2 ^= (!!(M_diag1 & shift)) * shift;
  M_diag1 ^= shift;
}

void DCHState::Sdag(unsigned target)
{
  //Add three to M_{ii}
  // 0 + 3 = 3; => 00 + 11 = 11;
  // 1 + 3 = 4 = 0; => 01 + 11 = 00;
  // 2 + 3 = 5 = 1; => 10 + 11 = 01;
  // 3 + 3 = 6 = 2; => 11 + 11 = 10;
  // Simply flip the ones bit, conditionally flip twos bit
  uint_t shift = (one << target);
  M_diag2 ^= ((!(M_diag1 & shift)) * shift);
  M_diag1 ^= (one << target);
}

void DCHState::Z(unsigned target)
{
  //Add two to M_{ii}
  M_diag2 ^= (one << target);
}

void DCHState::X(unsigned target)
{
  pauli_t p;
  p.X ^= (one << target);
  CommutePauli(p);
  s ^= p.X;
  omega.e += (hamming_parity(p.Z&s)*4);
  omega.e += (2*(p.e));
  omega.e = omega.e % 8;
}

void DCHState::Y(unsigned target)
{
  X(target);
  Z(target);
  omega.e = (omega.e + 6) % 8;
}

void DCHState::H(unsigned target)
{
  // Represent H_{a} = \frac{1}{2}\left(X_{a} + Z_{a}\right)
  pauli_t p, q;
  p.X ^= (one << target);
  CommutePauli(p);
  // Optimised commutation for the single Z pauli
  q.Z ^= G[target];
  q.X = (q.Z & v);
  q.Z = (q.Z & (~v));
  uint_t t = s^p.X;
  unsigned t_phase = (hamming_parity(t&p.Z) * 2U);
  t_phase += p.e;
  t_phase = t_phase % 4;
  uint_t u = s^q.X;
  unsigned u_phase = ( hamming_parity(s&q.Z) * 2U);
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
    b = (4-t_phase+2)%4;
  }
  // std::cout << "Global phase is: " << t_phase << std::endl;
  omega.e = (omega.e + 2*t_phase) %8;
  b = b%4;
  if(t==u)
  {
    // std::cout << "Found the same string" << std::endl;
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

void DCHState::CZ(unsigned control, unsigned target)
{
  if(control == target)
  {
    throw std::logic_error("Controlled operation cannot target the same qubit.");
  }
  M[target] ^= (one << control);
  M[control] ^= (one << target);
}

void DCHState::CX(unsigned control, unsigned target)
{
  if(control == target)
  {
    throw std::logic_error("Controlled operation cannot target the same qubit.");
  }
  isReadyGT = false;
  uint_t target_col = M[target];
  uint_t control_shift = (one << control);
  uint_t target_shift = (one << target);
  //Update the control col of M
  M[control] ^= target_col;
  bool target_bit = !!(M_diag1 & target_shift);
  for (unsigned i=0; i<n; i++)
  {
    // Update the control row of M
    // We use the target col as M is symmetric
    M[i] ^= (((target_col >> i) & one) * control_shift);
    //Update the control row of G^-1
    G_inverse[i] ^= ( !!(G_inverse[i] & target_shift)  * control_shift);
  }
  //Update control row/column with M_{t,t}
  M[control] ^= (target_bit * target_shift);
  M[target] ^= (target_bit * control_shift);
  // M_{c,c} += 2* M_{c,t}
  if ((target_col >> control) & one)
  {
    M_diag2 ^= control_shift;
  }
  // M_{c,c} += M_{t,t};
  if(!!(M_diag1 & control_shift) & target_bit)
  {
    M_diag2 ^= control_shift;
  }
  M_diag1 ^= (target_bit*control_shift);
  M_diag2 ^= ((!!(M_diag2 & target_shift))*control_shift);
  // Update the target column of G
  G[target] ^= G[control];
}

scalar_t DCHState::Amplitude(uint_t x)
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
  amp.e = (2*p.e);
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
  amp *= omega;
  return amp;
}

void DCHState::MeasurePauli(pauli_t P)
{
  CommutePauli(P);
  unsigned b = P.e;
  uint_t u = s^P.X;
  b += (hamming_parity(s&P.Z) * 2);
  b = b%4;
  UpdateSVector(s, u, b);
}

void DCHState::MeasurePauliProjector(std::vector<pauli_t>& generators)
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

void DCHState::TransposeG()
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

void DCHState::CommutePauliDC(pauli_t &p)
{
  unsigned phase = 0;
  if(p.X)
  {
    phase = (hamming_weight(p.X&M_diag1)%4) + (2*hamming_parity(p.X&M_diag2));
    uint_t offdiag_terms = zer;
    for (unsigned i=0; i<n; i++)
    {
      p.Z ^= ((one << i) * hamming_parity(p.X&M[i]));
      if((p.X>>i)& one)
      {
        offdiag_terms += hamming_weight(p.X&M[i]);
      }
    }
    //Compensate for the zero diagonal in M
    p.Z ^= (p.X & M_diag1);
    phase += ((!!(offdiag_terms & (one << 1))) * 2U);
    phase = (phase % 4);
  }
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
    x_temp ^= (shift * hamming_parity(p.X & G_inverse[i]));
  }
  p.X = x_temp;
  p.Z = z_temp;
}

void DCHState::CommutePauli(pauli_t& p)
{
  CommutePauliDC(p);
  uint_t x_temp = p.X;
  uint_t z_temp = p.Z;
  p.e = (p.e + (hamming_parity(p.X & p.Z & v)) * 2)%4;
  //Commute through H
  p.X = (z_temp & v) ^ (x_temp & (~v));
  p.Z = (x_temp & v) ^ (z_temp & (~v));
}

void DCHState::UpdateSVector(uint_t t, uint_t u, unsigned b)
{
  // std::cout << "Howdy from inside updateSvector!" << std::endl;
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
  // std::cout << "T string is :";
  // Print(t, n);
  // std::cout << std::endl;
  // std::cout << "U string is :";
  // Print(u, n);
  // std::cout << std::endl;
  uint_t nu0 = (t^u) & (~v);
  uint_t nu1 = (t^u) & v;
  unsigned q = 0;
  uint_t q_shift = zer;
  bool nu0Empty = true;
  // std::cout << "nu0: ";
  // Print(nu0, n);
  // std::cout << std::endl;
  // std::cout << "nu1: ";
  // Print(nu1, n);
  // std::cout << std::endl;
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
    nu0 &= (~q_shift);
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
    nu1 &= (~q_shift);
  }
  //Pick the string with the q-th bit equal to zero. If necessary, swap the strings.
  if (!!(t & q_shift))
  {
    // std::cout << "Swapping strings based on the qth bit" << std::endl;
    s=u;
    omega.e = (omega.e + 2*b) % 8;
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
  // std::cout << "H/S params: e1 "<< e1 << " e2:" << e2 << " e3:" << e3 << " e4:" << e4 << std::endl;
  //Split up based on if nu0 is empty
  if (nu0Empty)
  {
    // std::cout << "nu0 is empty" << std::endl;
    if (nu1)//if T2 != Identity
    {
      isReadyGT = false;
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
    if(e2) // Commute the s gate through the (updated) UC 
    {
      uint_t y = G_inverse[q];
      M_diag2 ^= (y & M_diag1);
      M_diag1 ^= y;
      for(unsigned i=0; i<n; i++)
      {
        M[i] ^= ((y>>i)&one)*y;
        M[i] &= ~(one << i); //Zero out the diagonal
      }
    }
  }
  else
  {
    // std::cout << "nu0 not empty" << std::endl;
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
      isReadyGT=false;
      for(unsigned i=0; i<n; i++)
      {
        if((nu0 >> i) & one)
        {
          //xor col q of G^{-1} into col i
          G_inverse[i] ^= G_inverse[q];
          uint_t shift = (one << i);
          for(unsigned j=0; j<n; j++)
          {
            //xor row i of G into row q
            G[j] ^= ((!!(G[j] & shift)) * q_shift);
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
  // std::cout << "V after:";
  // Print(v, n);
  // std::cout << std::endl;

  // update the scalar factor omega
  omega.e=(omega.e  + e1) % 8;
}

scalar_t DCHState::InnerProduct(uint_t &A_diag1, uint_t &A_diag2, std::vector<uint_t> &A)
{
    if(!isReadyGT)
    {
        TransposeG();
    }
    std::vector<uint_t> placeholder(n, zer);
    std::vector<uint_t> B(n, zer);
    uint_t B_diag1 = zer, B_diag2 = zer;
    std::vector<uint_t> J(n, zer);
    uint_t J_diag1 = zer, J_diag2 = zer;
    // Setup the matrix M-A, store it in B
    // (M-A) == (M + A*) => Flip 1s of A to 3s, then add...
    J_diag2 = (A_diag1 ^ A_diag2 ^ M_diag2 ^ (A_diag1 & M_diag1));
    J_diag1 = (A_diag1 ^ M_diag1);
    for(unsigned i=0; i<n; i++)
    {
        J[i] = M[i]^A[i]; //M is symmetric, use a column.
        J[i] ^= (J_diag1 & (one << i)); //Setup the ones on the diagonal
    }
    //Compute (M-A)G^{T}
    for(unsigned i=0; i<n; i++)
    {
      uint_t J_row = J[i]; //(M-A) is symmetric so grab a column
      uint_t shift = (one << i);
      for(unsigned j=0; j<n; j++)
      {
        unsigned weight = hamming_weight(J_row&GT[j]);
        placeholder[j] ^= ((weight & 1U) * shift);
      }
    }
    // Compute B:= G*((M-A)G^{T})
    for(unsigned i=0; i<n; i++)
    {
      uint_t G_row = GT[i]; //(M-A) is symmetric so grab a column
      uint_t shift = (one << i);
      for(unsigned j=0; j<n; j++)
      {
        bool parity = hamming_parity(G_row&placeholder[j]);
        B[j] ^= (parity * shift);
      }
      for(unsigned k=0; k<n; k++)
      {
        if((G[k]>>i)&one)
        {
          uint_t one_bit = (J_diag1 >> k) & one;
          if ((B_diag1 >> i) & one_bit)
          {
              B_diag2 ^= shift;
          }
          B_diag1 ^= one_bit * shift;
          B_diag2 ^= ((J_diag2 >> k) & one) * shift;
          for (size_t l=k+1; l<n; l++)
          {
            if ((J[l] >> k) & (G[l] >> i) & one)
            {
                B_diag2 ^= shift;
            }
          }
        }
      }
    }
    //Setup the quadratic form to evaluate the exponential sum
    QuadraticForm q(hamming_weight(v));
    q.Q = 4*(hamming_parity(s&v));
    unsigned col = 0;
    for(unsigned i=0; i<n; i++)
    {
      uint_t i_shift = (one << i);
      if(!!(s & i_shift))
      {
        q.Q = (q.Q + (!!(B_diag1 & i_shift))*2)%8;
        q.Q ^= (!!(B_diag2 & i_shift))*4;
        for(unsigned j=i+1; j<n; j++)
        {
          if((s >> j) & (B[j] >> i) & one)
          {
            q.Q ^= 4;
          }
        }
      }
      if (!!(v&i_shift))
      {
        uint_t col_shift = (one << col);
        q.D1 ^= (!!(B_diag1 & i_shift)) * col_shift;
        q.D2 ^= (!!((B_diag2 ^ s) & i_shift) ^ hamming_parity(B[i]&s)) * col_shift;
        unsigned row=0;
        for(unsigned j=0; j<n; j++)
        {
            if(!!(v&(one << j)))
            {
              q.J[col] |= (!!(B[i]&(one << j)))*(one << row);
              row++;
            }
        }
        col++;
      }
    }
    scalar_t amp = q.ExponentialSum();
    amp.p -= (n+q.n);
    amp *= omega;
    return amp;
}

scalar_t DCHState::InnerProduct(DCHState& other)
{
  //Compute the inner product between two DCH states
  // <phi_1|phi_2> = <s_1|H_1 C^{-1}_1 D^{-1}_1 D_2 C_2 H_2 |s_2>
  // Set J = D^{-1}_1D_2
  // Let J' = C^{-1}_1 J C_1
  // Set G = C^{-1} C_2
  // Rearrange <s_1| H_1 J' G H_2 |s_2>
  // Define |ip> = J' G H_2 |s_2>
  // Apply H_1 to |ip>
  // Compute <s_2 | ip' > and return
  if(this->n != other.n)
  {
    throw std::runtime_error("Cannot perform the inner product between states with "
                             "different numbers of qubits.");
  }
  DCHState *left, *right;
  if(hamming_weight(other.v) <= hamming_weight(this->v))
  {
    left = &other;
    right = this;
  }
  else
  {
    left = this;
    right = &other;
  }
  DCHState ip_state(right->n);
  ip_state.v = right->v;
  ip_state.s = right->s;
  ip_state.omega = right->omega;
  if(!(left->isReadyGT))
  {
    left->TransposeG();
  }
  //Compute the combined phase matrix
  uint_t diag_1 = (right->M_diag1 ^ left->M_diag1);
  uint_t diag_2 = (left->M_diag2 ^ left->M_diag1) ^ right->M_diag2 ^ (right->M_diag1 & left->M_diag1);
  std::vector<uint_t> placeholder(n, zer);
  std::vector<uint_t> J(n, zer);
  for(unsigned i=0; i<n; i++)
  {
    J[i] = (left->M[i] ^ right->M[i]);
    uint_t shift = (one << i);
    if(!!(diag_1 & shift))
    {
      J[i] ^= shift;
    }
  }
  //Commute the left CNOT matrix past the phase block
  for(unsigned i=0; i<n; i++)
  {
    uint_t J_row = J[i]; //(M-A) is symmetric so grab a column
    uint_t shift = (one << i);
    for(unsigned j=0; j<n; j++)
    {
      unsigned weight = hamming_weight(J_row&left->GT[j]);
      placeholder[j] ^= ((weight & 1U) * shift);
    }
  }
  // Compute B:= G*((M-A)G^{T})
  for(unsigned i=0; i<n; i++)
  {
    uint_t G_row = left->GT[i]; //(M-A) is symmetric so grab a column
    uint_t shift = (one << i);
    for(unsigned j=0; j<n; j++)
    {
      bool parity = hamming_parity(G_row&placeholder[j]);
      ip_state.M[j] ^= (parity * shift);
    }
    //Make sure the M matrix has zero diagonal...
    ip_state.M[i] &= ~(shift);
    for(unsigned k=0; k<n; k++)
    {
      if((left->G[k]>>i)&one)
      {
        uint_t one_bit = (diag_1 >> k) & one;
        if ((ip_state.M_diag1 >> i) & one_bit)
        {
            ip_state.M_diag2 ^= shift;
        }
        ip_state.M_diag1 ^= one_bit * shift;
        ip_state.M_diag2 ^= ((diag_2 >> k) & one) * shift;
        for (size_t l=k+1; l<n; l++)
        {
          if ((J[l] >> k) & (left->G[l] >> i) & one)
          {
              ip_state.M_diag2 ^= shift;
          }
        }
      }
    }
  }
  //Compute the combined G matrix
  for(unsigned i=0; i<n; i++)
  {
    uint_t shift = (one << i);
    for(unsigned j=0; j<n; j++)
    {
      ip_state.G[j] ^= hamming_parity(right->GT[i]&left->G_inverse[j])*shift;
      ip_state.G_inverse[j] ^= hamming_parity(left->GT[i]&right->G_inverse[j])*shift;
    }
  }
  // Commute Hadamards through to the IP state
  for(unsigned i=0; i<n; i++)
  {
    if((left->v>> i) & one)
    {
      ip_state.H(i);
    }
  }
  // Get the right global phase for the left hand state
  scalar_t left_amp(left->omega);
  left_amp.conjugate();
  // Calculate the amplitude, correct and return.
  scalar_t amp = ip_state.Amplitude(left->s);
  amp *= left_amp;
  return amp;
}

}
#endif