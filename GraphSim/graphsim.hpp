// graphsim.h

/*!\mainpage
This is 'graphsim', a simulator for stabilizer quantum cuircuits using the
graph state formalism.

graphsim version 0.10, dated 2005-Jan-27

(c) Simon Anders (sanders@fs.tum.de), University of Innsbruck

See the article

  S. Anders, H. J. Briegel: 
  Fast simulation of Stabilizer Circuits using a Graph States Formalism
  quant-ph/0504117

for a description of this software.

----

// This file is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.

// This file is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this file; see the file COPYING.  If not, browse to 
// http://www.fsf.org/licenses/gpl.html

*/

/*!\file
This header file defines the main interface of graphsim.

(c) Simon Anders, University of Innsbruck, 2005
released under GPL.

Note: If you have trouble compiling this, please note:
This file uses the "hash_set" template, which is an extension
to the Standard C++ Library and the Standard Template Library,
specified by SGI. Most C++ compilers have this template.
For GNU C++, the header file <ext/hash_set> hasa to be included
and hash_set has to be prefixed with namespace __gnu_cxx.
If you use another compiler, you might have to change the include
file and the namespace identifier.
*/

#ifndef GRAPHSIM_H
#define GRAPHSIM_H

#include <cassert>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <stdexcept>
#include <unordered_set>
#include <utility>
#include <ctime>


#include "loccliff.hpp"

#include <cstring>
#include <cstdlib>

namespace GraphSim
{
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

/*! All vertices in a graph state are numbered beginning with 0. To specify
auch an index, the type VertexIndex (which is just unsigned long) is always
used */
typedef unsigned long VertexIndex;

/*! A GraphRegister object maintains a list of its vertices (qubits), each 
described by an object of the class QubitVertex described here.*/
struct QubitVertex {
   /*!byprod is the vertex operator (VOp) associated with the qubit (the name 
   stems from the term 'byproduct operator' used for the similar concept in 
   the one-way quantum computer.*/
   LocCliffOp byprod;
   /*! neigbors is the adjacency list for this vertex */
   std::unordered_set<VertexIndex> neighbors;
   /*! Upon construction, a qubit vertex is initialised with the Hadamard
   operation as VOp, and with wmpty neighbor list. This makes it represent
   a |0>. */
   QubitVertex (void)
     : byprod (lco_H) {};
};
   
/*! This structure is only for internal use for the cphase functions. */
struct ConnectionInfo {
   bool wasEdge;
   bool non1;
   bool non2;
};

//! A quantum register.
/*! GraphRegister is the central class of graphsim. It represents a register of qubits
that can be entangled with each other. It offers functions to initialize the register,
let gates operate on the qubits, do measurements and print out the state. */
class GraphRegister {
  public:
   /*! This vector stores all the qubits, represented as QubitVertex objects. The index
   of the vector is usually taken as of type VertexIndex. */
   std::vector<QubitVertex> vertices;
   GraphRegister (VertexIndex numQubits, int randomize = -1);
   GraphRegister (GraphRegister& gr);
   ~GraphRegister () {};
   void local_op (VertexIndex v, LocCliffOp o);
   void hadamard (VertexIndex v);
   void phaserot (VertexIndex v);
   void bitflip (VertexIndex v);
   void phaseflip (VertexIndex v);
   void cphase (VertexIndex v1, VertexIndex v2);      
   void cnot (VertexIndex vc, VertexIndex vt);
   int measure (VertexIndex v, LocCliffOp basis = lco_Z, 
      bool* determined = NULL, int force = -1);
   Stabilizer & get_full_stabilizer (void) const;
   void invert_neighborhood (VertexIndex v);
   void print_adj_list (std::ostream& os = std::cout) const;
   void print_adj_list_line (std::ostream& os, VertexIndex i) const;
   void print_stabilizer (std::ostream& os = std::cout) const;
  private:
   void add_edge (VertexIndex v1, VertexIndex v2);
   void del_edge (VertexIndex v1, VertexIndex v2);
   void toggle_edge (VertexIndex v1, VertexIndex v2);
   int graph_Z_measure (VertexIndex v, int force = -1);
   int graph_Y_measure (VertexIndex v, int force = -1);
   int graph_X_measure (VertexIndex v, bool* determined = NULL, int force = -1);
   void toggle_edges (const std::unordered_set<VertexIndex> vs1, 
      const std::unordered_set<VertexIndex> vs2);      
   bool remove_byprod_op (VertexIndex v, VertexIndex use_not);
   void cphase_with_table (VertexIndex v1, VertexIndex v2);
   ConnectionInfo getConnectionInfo (VertexIndex v1, VertexIndex v2);  
};

/*! As we often iterate over sublists of GraphRegister::vertices, this
iterator typedef is a handy abbreviation. */
typedef std::vector<QubitVertex>::iterator VertexIter;
/*! Another iterator, this one for the adjacency lists QubitVertex::neigbors,
and subsets. */
typedef std::unordered_set<VertexIndex>::iterator VtxIdxIter;

/*! A constant version of VertexIter */
typedef std::vector<QubitVertex>::const_iterator VertexIterConst;
/*! A constant version of VtxIdxIter */
typedef std::unordered_set<VertexIndex>::const_iterator VtxIdxIterConst;

/*! Apply the local (i.e. single-qubit) operation o on vertex v. */
inline void GraphRegister::local_op (VertexIndex v, LocCliffOp o) {
   vertices[v].byprod = o * vertices[v].byprod;
}

/*! Apply a Hadamard gate on vertex v */
inline void GraphRegister::hadamard (VertexIndex v) {
   local_op (v, lco_H);
} 

/*! Apply a phaserot gate on vertex v. Phaserot means the gate
S = |0><0| + i |1><1|. */
inline void GraphRegister::phaserot (VertexIndex v) {
   local_op (v, lco_S);
} 

/*! Apply a bitflip gate (i.e. a Pauli X) on vertex v */
inline void GraphRegister::bitflip (VertexIndex v) {
   local_op (v, lco_X);
} 

/*! Apply a phaseflip gate (i.e. a Pauli Z) on vertex v */
inline void GraphRegister::phaseflip (VertexIndex v) {
   local_op (v, lco_Z);
} 

//!Random coin toss. Change here to insert your favorite RNG.
int bool_rand (void) {
  return random () > RAND_MAX/2;
}

//! Instantiate a quantum register with 'numQubits' qubits, initally all in state |0>. 
/*! If randomize > -1 the RNG will be seeded with the current time plus the value 
of randomize. (Otherwise, it is not seeded.) That the value of randomize is 
added to the seed is useful in parallel processing settings where you want to
ensure different seeds. (If you call this from Python, remember, that Python's 
RNG is not seeded.) */
GraphRegister::GraphRegister (VertexIndex numQubits, int randomize)
   : vertices (numQubits) 
{
   if (randomize > -1) {
      srandom (time(NULL) + randomize);
   }
}

//! Copy constructor
/*! Clones a register */
GraphRegister::GraphRegister (GraphRegister& gr)
   : vertices (gr.vertices)
{
}

//! Add an edge to the graph underlying the state.
void GraphRegister::add_edge (VertexIndex v1, VertexIndex v2) {
   assert (v1 != v2);
   vertices[v1].neighbors.insert (v2);
   vertices[v2].neighbors.insert (v1);
   //D cerr << "adding edge " << v1 << " - " << v2 << endl;
}

//! Delete an edge to the graph underlying the state.
void GraphRegister::del_edge (VertexIndex v1, VertexIndex v2) {
   vertices[v1].neighbors.erase (v2);
   vertices[v2].neighbors.erase (v1);
}

//! Toggle an edge to the graph underlying the state.
/*! (i.e. add it if not present, and delete it if present.) */
void GraphRegister::toggle_edge (VertexIndex v1, VertexIndex v2) 
{
   int n1 = vertices[v1].neighbors.erase (v2);
   if (n1 == 1) {
      vertices[v2].neighbors.erase (v1);
   } else {
      assert (n1 == 0);
      add_edge (v1, v2);
   }
}

//! Create the Stabilizer of the state. 
/*! This is useful to print out the stabilizer (or to compare with CHP). 
You can also use print_stabilizer.*/
Stabilizer& GraphRegister::get_full_stabilizer (void) const
{
   std::unordered_set<VertexIndex> all_qubits; 
   for (VertexIterConst i = vertices.begin(); i != vertices.end(); i++) {
      all_qubits.insert (i-vertices.begin());
   }
   Stabilizer *s = new Stabilizer (*this, all_qubits);
   return *s;
}

//! Prints out the description of the current state 
/*! in terms of adjacency lists of the graph and the VOps.*/
void GraphRegister::print_adj_list (std::ostream& os) const
{
   for (VertexIndex i = 0; i < vertices.size(); i++) {
      print_adj_list_line (os, i);
   }  
}

//! Prints the line for Vertex i in the adjacency list representation of the state.
void GraphRegister::print_adj_list_line (std::ostream& os, VertexIndex i) const 
{
   os << "Vertex " << i << ": VOp "
      <<  vertices[i].byprod.get_name() << ", neighbors ";
   for (VtxIdxIterConst j = vertices[i].neighbors.begin(); 
         j != vertices[i].neighbors.end(); j++) {
      os << *j << " ";
   }
   os << std::endl;
}

//! Print the current state in stabilizer representation.
void GraphRegister::print_stabilizer (std::ostream &os) const
{
   get_full_stabilizer().print (os);
}

/*! Use the cphase look-up table. This is called by cphase after VOps that do not
commute with the cphase gate have been removed as far as possible. */
void GraphRegister::cphase_with_table (VertexIndex v1, VertexIndex v2)
{
   // Structure of the table:
   // first index: whether there was an edge between the operands before
   //   (0=no, 1= yes)
   // second and third index: byprod op of v1 and v2
   // third index: information to obtain:
   //    0= whether after the cphase there is an edges
   //    1,2= new values of the byprod ops of v1 and v2Id
   static const unsigned short cphase_tbl[2][24][24][3] =
      #include "cphase.tbl"
   ;   
   ConnectionInfo ci = getConnectionInfo (v1, v2);
   unsigned op1 = vertices[v1].byprod.op;
   unsigned op2 = vertices[v2].byprod.op;
   // The table must only be used if a vertex has either no
   // non-operand neighbors, or a diagonal byprod op
   assert ((!ci.non1) || vertices[v1].byprod.is_diagonal());
   assert ((!ci.non2) || vertices[v2].byprod.is_diagonal());
   if (cphase_tbl[ci.wasEdge][op1][op2][0]) {
      add_edge (v1, v2);
   } else {
      del_edge (v1, v2);
   }
   vertices[v1].byprod.op = cphase_tbl[ci.wasEdge][op1][op2][1];
   vertices[v2].byprod.op = cphase_tbl[ci.wasEdge][op1][op2][2];
   // The condition above must also hold afterwards:  
   ci = getConnectionInfo (v1, v2);
   assert ((!ci.non1) || vertices[v1].byprod.is_diagonal());
   assert ((!ci.non2) || vertices[v2].byprod.is_diagonal());
}   

/*! Check whether the qubits are connected to each other and to non-operand vertices.*/
ConnectionInfo GraphRegister::getConnectionInfo 
      (VertexIndex v1, VertexIndex v2) {
   ConnectionInfo ci;
   ci.wasEdge = 
     vertices[v1].neighbors.find(v2) != vertices[v1].neighbors.end();
   if (! ci.wasEdge) {
      ci.non1 = vertices[v1].neighbors.size() >= 1; 
      ci.non2 = vertices[v2].neighbors.size() >= 1;
   } else {
      ci.non1 = vertices[v1].neighbors.size() >= 2; 
      ci.non2 = vertices[v2].neighbors.size() >= 2;
   }
   return ci;
}
   
//! Do a conditional phase gate between the two qubits.
void GraphRegister::cphase (VertexIndex v1, VertexIndex v2) 
{
   // If there are non-operand neighbors, we can use neighborhood inversion
   // to remove the byprod operators. 
   // These will store whether the operand vertices have nonoperand neighbors.
   ConnectionInfo ci = getConnectionInfo (v1, v2);
   
   if (ci.non1) {
      remove_byprod_op (v1, v2);
   }
   ci = getConnectionInfo (v1, v2);
   if (ci.non2) {
      remove_byprod_op (v2, v1);
   }
   ci = getConnectionInfo (v1, v2);
   if (ci.non1 && !vertices[v1].byprod.is_diagonal()) {
      // this can happen if v1 was first skipped
      remove_byprod_op (v1, v2);
   }
   cphase_with_table (v1, v2);
}   


//! Do a controlled not gate between the vertices vc (control) and vt (target).         
void GraphRegister::cnot (VertexIndex vc, VertexIndex vt) {
   hadamard (vt);
   cphase (vc, vt);
   hadamard (vt);
}

//! This structure is needed only by toggle_edges
struct Edge: public std::pair<VertexIndex, VertexIndex> {
   Edge (const VertexIndex a, const VertexIndex b) {
      if (a < b) {
         first = a; second = b;
      } else {
         first = b; second = a;
      }
   }
};

//! This structure is needed only by toggle_edges
struct edge_hash {
   size_t operator() (const Edge& e) const
   {
      return e.first << 16 ^ e.second;
   };
};

//! Toggles the edges between the vertex sets vs1 and vs2.
/*! The function takes extra care not to invert an edge twice. If vs1 and
vs2 are disjunct, this cannot happen and we do not need the function.
If vs1 == v2s, we can do without, too. */
void GraphRegister::toggle_edges (const std::unordered_set<VertexIndex> vs1, 
   const std::unordered_set<VertexIndex> vs2)
{
   std::unordered_set<Edge, edge_hash> procd_edges;
   for (VtxIdxIterConst i = vs1.begin(); i != vs1.end(); i++) {
      for (VtxIdxIterConst j = vs2.begin(); j != vs2.end(); j++) {
         if ((*i != *j) && 
               (procd_edges.find (Edge (*i, *j)) == procd_edges.end())) {
            procd_edges.insert (Edge (*i, *j));
            toggle_edge (*i, *j);
         }
      }
   }
}

//! Measure the bare graph state in the Z basis.
int GraphRegister::graph_Z_measure (VertexIndex v, int force) 
{
   int res;
   if (force == -1) {
      res = bool_rand ();
   } else {
      res = force;
   }
   std::unordered_set<VertexIndex> nbg = vertices[v].neighbors;
   for (VtxIdxIter i = nbg.begin(); i != nbg.end(); i++) {
      del_edge (v, *i);
      if (res) {
         vertices[*i].byprod = vertices[*i].byprod * lco_Z;
      }
   }
   if (! res) {
      vertices[v].byprod = vertices[v].byprod * lco_H;
   } else {
      vertices[v].byprod = vertices[v].byprod * lco_X * lco_H;
   }
   return res;   
}

//! Measure the bare graph state in the Y basis.
int GraphRegister::graph_Y_measure (VertexIndex v, int force)
{
   int res;
   if (force == -1) {
      res = bool_rand ();
   } else {
      res = force;
   }
   std::unordered_set<VertexIndex> vnbg = vertices[v].neighbors;
   for (VtxIdxIter i = vnbg.begin(); i != vnbg.end(); i++) {
      if (res) {
         vertices[*i].byprod = vertices[*i].byprod * lco_spiZ;
      } else {
         vertices[*i].byprod = vertices[*i].byprod * lco_smiZ;
      }
   }
   vnbg.insert (v); // Now, vnbg is the set of v and its neighbours.
   for (VtxIdxIter i = vnbg.begin(); i != vnbg.end(); i++) {
      for (VtxIdxIter j = i; j != vnbg.end(); j++) {
         if (i != j) {
           toggle_edge (*i, *j);
         }   
      }
   }
   if (! res) {
      // Messergebnis: +|0y>
      vertices[v].byprod = vertices[v].byprod * lco_S;
   } else {
      // Messergebnis: -|0y>
      vertices[v].byprod = vertices[v].byprod * lco_S.herm_adjoint();
   }
   return res;   
}

//! Measure the bare graph state in the X basis.
int GraphRegister::graph_X_measure (VertexIndex v, bool* determined, 
      int force) 
{
   if (vertices[v].neighbors.size () == 0) {
      //not entangled qubit => result always 0
      if (determined != NULL)
         *determined = true;
      return 0;
   }
   if (determined != NULL)
      *determined = false;
   // entangled qubit => let's get on with the complicated procedure
   // throw a die:
   int res;
   if (force == -1) {
      res = bool_rand ();
   } else {
      res = force;
   }
   VertexIndex vb = *vertices[v].neighbors.begin(); // the choosen vertex
   //D cerr << "vb = " << vb << endl;
   // preparation step: store the neighborhood of v and vb
   std::unordered_set<VertexIndex> vn = vertices[v].neighbors;
   std::unordered_set<VertexIndex> vbn = vertices[vb].neighbors;
   // First, put the byproduct ops: 
   if (! res) {
      // measured a |+>:
      // lco_spiY on vb
      vertices[vb].byprod = vertices[vb].byprod * lco_spiY;
      // Z on all in nbg(v) \ nbg(vb) \ {vb}
      for (VtxIdxIter i = vertices[v].neighbors.begin(); 
            i != vertices[v].neighbors.end(); i++) {
         if (*i != vb && 
               vertices[vb].neighbors.find(*i) == vertices[vb].neighbors.end()) {  
            vertices[*i].byprod = vertices[*i].byprod * lco_Z;
         }
      }
   } else {
      // measured a |->:
      // lco_smiY on vb, and lco_Z on v:
      vertices[vb].byprod = vertices[vb].byprod * lco_smiY;
      vertices[v].byprod = vertices[v].byprod * lco_Z;
      // Z on all in nbg(vb) \ nbg(v) \ {v}
      for (VtxIdxIter i = vertices[vb].neighbors.begin(); 
            i != vertices[vb].neighbors.end(); i++) {
         if (*i != v && 
               vertices[v].neighbors.find(*i) == vertices[v].neighbors.end()) {  
            vertices[*i].byprod = vertices[*i].byprod * lco_Z;
         }
      }
   }     
   // Toggling the edges in three steps
   // STEP 1: complement with Edges (nbg(v), nbg(vb)):
   toggle_edges (vn, vbn);
   // STEP 2: complement with the complete subgraph induced by the 
   // intersection of nbg(v) and nbg(vb):
   // First, make the intersection
   std::unordered_set<VertexIndex> isc;
   for (VtxIdxIter i = vn.begin(); i != vn.end(); i++) {
      if (vbn.find(*i) != vbn.end()) {
         isc.insert (*i);
      }
   }
   // Now, toggle the edges
   for (VtxIdxIter i = isc.begin(); i != isc.end(); i++) {
      for (VtxIdxIter j = i; j != isc.end(); j++) {
         if (*i != *j) {
            toggle_edge (*i, *j);
         }
      }
   }
   // STEP 3: Toggle all edges from vb to nbg(v) \ {vb}
   for (VtxIdxIter i = vn.begin(); i != vn.end(); i++) {
      if (*i != vb) {
         toggle_edge (vb, *i);
      }
   }
   return res;   
}

//! Measure qubit v in basis 'basis'.
/* For basis, pass a LocCliffOp object which has to be equal to one of lco_X, lco_Y
or lco_Z. If you want to now, whether the result was choosen at random or determined
by the state, pass a bool pointer in which this information will be written.
If you want to force the result to be a certain value, pass 0 or 1 to 'force'. This only
works, if the result is not determined. If it is, 'force' is ignored.*/
int GraphRegister::measure (VertexIndex v, LocCliffOp basis, bool* determined,
      int force) 
{
   assert (basis.op >= lco_X.op && basis.op <= lco_Z.op);
   if (determined != NULL) {
      *determined = false;
   }
   assert (force >= -1 && force <= 1);
   LocCliffOp basis_orig = basis;
   RightPhase rp = basis.conjugate (vertices[v].byprod.herm_adjoint());
   assert (rp == rp_p1 || rp == rp_m1);
   if (force != -1 && rp == rp_m1) {
      force = force ^ 0x01;
   }
   int res;
   switch (basis.op) {
      case 1 /* lco_X */: res = graph_X_measure (v, determined, force); break;
      case 2 /* lco_Y */: res = graph_Y_measure (v, force); break;
      case 3 /* lco_Z */: res = graph_Z_measure (v, force); break;
      default: exit (1);
   }
   if (rp == rp_m1) {
      res = ! res;
   } else {
      assert (rp == rp_p1);
   }
   // check: the measured vertex should be singled out:
   assert (vertices[v].neighbors.size() == 0);
   // Check that the vertex is now in the correct eigenstate:
   LocCliffOp assert_op = lco_X;
   assert (assert_op.conjugate (vertices[v].byprod) == (res ? rp_m1 : rp_p1));
   assert (assert_op == basis_orig);
   return res;
}
  
//! Do a neighborhood inversion (i.e. local complementation) about vertex v.
/* This changes the state's graph representation but not the state itself, as the
necessary correction to the VOps are applied. */
void GraphRegister::invert_neighborhood (VertexIndex v)
{
   // Invert the neighborhood:
   std::unordered_set<VertexIndex> vn = vertices[v].neighbors;
   for (VtxIdxIter i = vn.begin(); i != vn.end(); i++) {
      for (VtxIdxIter j = i; j != vn.end(); j++) {
         if (*i != *j) {
            //cerr << "toggling " << *i << "," << *j << endl;
            toggle_edge (*i, *j);
         }
      }
      // and adjust the local Cliffords:
      vertices[*i].byprod = vertices[*i].byprod * lco_spiZ.herm_adjoint();
   }
   // finally, adjust the local Clifford of v:
   vertices[v].byprod = vertices[v].byprod * lco_smiX.herm_adjoint();
}

//! Do neighborhood inversions to reduce the VOp of vertex v to the identity.
/* 'avoid' is avoided as swapping partner, i.e. the swapping partner will not 
be 'avoid' unless this is the only neighbor of v. If no neighbor is available, 
the program exits with error message.*/
bool GraphRegister::remove_byprod_op (VertexIndex v, VertexIndex avoid)
{
   // A lookup table on how any LC operator can be composed from them
   // generators lco_spiZ (denoted 'V') and lco_smiX (denoted 'U'):
   // (The lookup table was generated from sqrt-gen.nb and copied here manually.)
   const char * comp_tbl [24] =
     {"UUUU", "UU", "VVUU", "VV", 
      "VUU", "V", "VVV", "UUV",
      "UVU", "UVUUU", "UVVVU", "UUUVU",
      "UVV", "VVU", "UUU", "U",
      "VVVU", "UUVU", "VU", "VUUU", 
      "UUUV", "UVVV", "UV", "UVUU"};
   // Of course, we need a neighborhood
   if (vertices[v].neighbors.size() == 0) {
      throw std::runtime_error("remove_byprod_op: called with isolated vertex. Aborting.\n");      
   }
   // This will be the swapping partner:
   VertexIndex vb = * vertices[v].neighbors.begin ();
   if (vb == avoid) {
      // Is there an alternative to 'avoid'? If so, use it.
      VtxIdxIter vbi = vertices[v].neighbors.begin ();
      vbi++;
      if (vbi != vertices[v].neighbors.end()) {
         vb = *vbi;
      }
   }   
   const char * comp = comp_tbl[vertices[v].byprod.op];
   for (int pos = strlen(comp)-1; pos >= 0; pos--) {
      if (comp[pos] == 'U') {
         // A U will vanish if we do an inversion on v
         invert_neighborhood (v);
      } else {
         assert (comp[pos] == 'V');
         // For this we need to invert on a neighbor of v
         invert_neighborhood (vb);
      }
   }
   // Now, we should have lco_Id left
   assert (vertices[v].byprod == lco_Id); 
   return true;
}

Stabilizer::Stabilizer (const VertexIndex numQubits_):
   paulis (numQubits_, std::vector<LocCliffOp> (numQubits_, lco_Id)),
   rowsigns (numQubits_, rp_p1),
   vtxidx (numQubits_, 0)
{
   numQubits = numQubits_;
}


Stabilizer::Stabilizer (const GraphRegister& gr, 
      const std::unordered_set<VertexIndex>& qubits):
   paulis (qubits.size(), std::vector<LocCliffOp> (qubits.size(), lco_Id)),
   rowsigns (qubits.size()),
   vtxidx (qubits.size())
{
   numQubits = qubits.size ();
   // Build the graph adjacency matrix with Z's and X's in the diagonal
   // and apply the local Clifford unitaries:
   int in = 0;
   for (VtxIdxIter i = qubits.begin(); i != qubits.end(); i++, in++) {
      rowsigns[in] = RightPhase (0);
     vtxidx[in] = *i;
      int jn = 0;
      for (VtxIdxIter j = qubits.begin(); j != qubits.end(); j++, jn++) {
         if (i==j) {
            paulis[in][jn] = lco_X;
         } else {
           if (gr.vertices[*i].neighbors.find(*j) !=
                 gr.vertices[*i].neighbors.end()) {
              paulis[in][jn] = lco_Z;
           } else {
              paulis[in][jn] = lco_Id;
           }
         }
         // Now the local Clifford unitaries:
         conjugate (in, jn, gr.vertices[*j].byprod);
      }   
   }
}

Stabilizer::Stabilizer (struct QState * qs) :
   paulis (qs->n, std::vector<LocCliffOp> (qs->n, lco_Id)),
   rowsigns (qs->n),
   vtxidx (qs->n)
{
   static const unsigned char optbl[4] = {0, 1, 3, 2};
   numQubits = qs->n;
   for (int i = 0; i < (int)numQubits; i++) {
      rowsigns[i] = RightPhase (qs->r[numQubits + i]);
     vtxidx[i] = i;
      for (int j = 0; j < (int)numQubits; j++) {
         bool xhere = ((qs->x [numQubits+i] [j >> 5]) & (1 << (j & 0x1f))) > 0;
         bool zhere = ((qs->z [numQubits+i] [j >> 5]) & (1 << (j & 0x1f))) > 0;
         paulis[i][j] = LocCliffOp (optbl [(zhere<<1) | xhere]);
      }
   }
}

void Stabilizer::add_row (unsigned target, unsigned addend)
{
   //D cerr << "adding row " << addend << " to row " << target << endl;
   for (unsigned col = 0; col < numQubits; col++) {
      rowsigns[target] = rowsigns[target] + 
        LocCliffOp::mult_phase (paulis[target][col], paulis[addend][col]);
      paulis[target][col] = paulis[target][col] * paulis[addend][col];
   }
}

void Stabilizer::conjugate (unsigned row, unsigned col, 
   const LocCliffOp trans)
{
  rowsigns[row] = rowsigns[row] + paulis[row][col].conjugate (trans);   
}

void Stabilizer::conjugate_column (unsigned col, const LocCliffOp trans)
{
   for (unsigned row = 0; row < numQubits; row++) {
      conjugate (row, col, trans);
   }
}

void Stabilizer::print (std::ostream &os) const
{
   for (unsigned i = 0; i < numQubits; i++) {
      os << rowsigns[i].get_name() << " ";
      for (unsigned j = 0; j < numQubits; j++) {
         os << paulis[i][j].get_name().substr(0,1) << " ";
      }
      os << std::endl;
   }
}

}
#endif //GRAPHSIM_H
