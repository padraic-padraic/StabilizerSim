// stabilizer.h

/*!\file 
defines the Stabilizer class */

#ifndef STABILIZER_H
#define STABILIZER_H

#include <iostream>
#include <vector>
#include <unordered_set>
#include <cassert>
#include "loccliff.hpp"

namespace GraphSim
{
// The following struct declaration is quoted in verbatim from Scott
// Aaronson's CHP program. (Scott, I hope, you don't mind.) It is used
// here to implement the functionality of comparing a stabilizer state
// in CHP with a state in graphsim. This is used to compare the action
// of the two programs to ensure correctness.
struct QState
{
   long n;         // # of qubits
   unsigned long **x; // (2n+1)*n matrix for stabilizer/destabilizer x bits (there's one "scratch row" at
   unsigned long **z; // (2n+1)*n matrix for z bits                                                 the bottom)
   int *r;         // Phase bits: 0 for +1, 1 for i, 2 for -1, 3 for -i.  Normally either 0 or 2.
   unsigned long pw[32]; // pw[i] = 2^i
   long over32; // floor(n/8)+1
}; 

//forward declarations, definition in graphsim.h:

typedef unsigned long VertexIndex;
struct QubitVertex;
class GraphRegister;

//! A stabilizer describing a quantum state.
/* This class is used at the moment only by GraphRegister::print_stabilizer and for
the compare method (when comparing with CHP). */
struct Stabilizer {
   VertexIndex numQubits;
   std::vector<std::vector<LocCliffOp> > paulis;
   std::vector<RightPhase> rowsigns;
   std::vector<VertexIndex> vtxidx;
   Stabilizer (const VertexIndex numQubits_);
   Stabilizer (const GraphRegister& gr, const std::unordered_set<VertexIndex>& qubits);
   Stabilizer (QState * qs);
   void add_row (unsigned target, unsigned addend);
   void conjugate (unsigned row, unsigned col, const LocCliffOp trans);
   void conjugate_column (unsigned col, const LocCliffOp trans);
   void print (std::ostream &os = std::cout) const;
};

} //end namespace GraphSim

#endif //STABILIZER_H
