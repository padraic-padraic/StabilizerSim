#include "stabilizer.hpp"
#include "graphsim.hpp"

namespace GraphSim
{
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
