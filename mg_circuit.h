// ----------------------------------------------------------------//
// Filename : mg_ckt.h
// Author : Ting Yu <tingyu1@illinois.edu>
//
// declaration of multi_grid Circuit class
// used to construct the circuit network
// ----------------------------------------------------------------//
// - Ting Yu - Tue Feb 8 5:45 pm 2011
//   * added the ostream<< func in .h file
// - Zigang Xiao - Tue Jan 18 21:46:13 CST 2011
//   * added solve() and related function
// - Zigang Xiao - Sun Jan 16 16:04:03 CST 2011
//   * added this log

#ifndef __MG_CIRCUIT_H__
#define __MG_CIRCUIT_H__

#include <string>
#include <vector>
#include <tr1/unordered_map>
#include <cmath>
#include "circuit.h"
using namespace std;

class MG_Circuit{
public:
	MG_Circuit();
	~MG_Circuit();

	int LEVEL; // multigrid levels
	vector<Circuit*> mg_ckt; // define circuit vector
	// add unordered map for VDD_set between nbr layers
	
	// functions
	// build up mg_ckts
	void build_mg_ckt(Circuit *ckt, int layer);
	Circuit * build_one_layer_circuit(Circuit *ckt);
	void set_nbr_nets(Node *nd, Node *&nd_c, Circuit *ckt,
		Circuit *&coarse_ckt);
	// solve the circuit from coarse to fine
	void solve_mg_ckt(Circuit *ckt);
};

#endif

