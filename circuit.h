// ----------------------------------------------------------------//
// Filename : circuit.h
// Author : Xiao Zigang <zxiao2@illinois.edu>
//
// declaration of Circuit class
// used to construct the circuit network
// ----------------------------------------------------------------//
// - Ting Yu - Tue Feb 8 5:45 pm 2011
//   * added the ostream<< func in .h file
// - Zigang Xiao - Tue Jan 18 21:46:13 CST 2011
//   * added solve() and related function
// - Zigang Xiao - Sun Jan 16 16:04:03 CST 2011
//   * added this log

#ifndef __CIRCUIT_H__
#define __CIRCUIT_H__

#include <string>
#include <vector>
#include <queue>
#include <tr1/unordered_map>
#include <map>
#include <list>
#include <cmath>
#include "cholmod.h"
#include "global.h"
#include "node.h"
#include "net.h"
#include "vec.h"
#include "triplet.h"
#include "block.h"
#include <set>
#include "circularqueue.h"
using namespace std;
using namespace std::tr1;

typedef vector<double> DoubleVector;
typedef vector<Net *> NetPtrVector;
typedef list<Net *> NetPtrList;
typedef vector<Node *> NodePtrVector;
typedef NetPtrVector NetList;

struct VoltageLessThan{
	inline bool operator() (const Node* node1, const Node* node2) {
    	return node1->get_value() < node2->get_value();
	}
};
// functor of translating Node * to void *
namespace std{ 
	namespace tr1{
		template<> struct hash< Node * >{
			size_t operator()( const Node * x ) const {
				return hash< const char* >()( (char*) x );
			}
		};
	}
}

class Circuit{
public:
	Circuit(string name="");
	~Circuit();
	void check_sys() const;
	friend class Block;
	// can be written as inline to speed up
	Node * get_node(string name);
	Net * get_net(string name);
	string get_name() const;

	static size_t get_total_num_layer();

	// add a node into nodelist
	bool add_node(Node * nd);

	// add a net into netset
	bool add_net(Net * net);
	bool add_net_pad(Net * net);

	bool has_node(string name) const;
	bool has_net(string name) const;

	// sort nodes according to predefined order
	void sort_nodes();

	// solve for node voltage
	void solve(double Frozen_T);
	void solve_coarse(double Frozen_T);
	// calculate power consumption for each node
	void compute_power();
	
	void set_blocklist(Node * nd);

	static void set_parameters(double, double, double, size_t, int);
	static void get_parameters(double&, double&, double&, size_t&, int&);

	friend ostream & operator << (ostream & os, const Circuit & ckt);
	friend class Parser;

	// C style output
	void print();

	/////////////////// function for Pad relocation //////
	//////// variables ////////
	double max_IRdrop;
	double th_IRdrop;
	// step is used to divide max_IRdrop into several classes
	// min_IRdrop should be 0
	double step;
	vector <Node*> VDD_set;
	vector <Node*> VDD_candi_set;
	set<Node*, VoltageLessThan> CriticalNodes;
	/////////// functions //////////
	void print_power();
	double optimize_pad_assign();
	double optimize_pad_assign_new();
	void clear_candi_visit();
	void rebuild_voltage_nets();
	size_t rebuild_voltage_nets_one_move(Node *rm_node, 
		Node *add_node, Net *&rm_net, Net *&add_net);

	// output the format for matlab, to plot the 
	// distribution figure;
	void print_matlab();
	void pad_set_init();

	// functions for new optimization 
	// method
	Node* find_min_IRdrop_pad(size_t & min_index);
	Node* find_max_IRdrop_candi();
	double area_IRdrop(Node *nd);

	//locate the max IR drop numbers
	double locate_maxIRdrop();
	void locate_thIRdrop();
	void random_init_iter();
	void RANSAC_init();
	void opti_SA();
	void Mean_shift_move();
	void Mean_shift_one_move(Node *center, 
	  size_t center_index, Node *&new_center, 
	  size_t &new_candi_index);
	Node* find_max_IRdrop_candi(Node *rm_pad, size_t Radius);
	void build_criticalNodes();
	double SACost();
	bool acceptProb(double p);
	double penalty(double v, double vworst);
	void stamp_rhs_SA(double* b);
	void solve_GS(double *b);
	void solve_GS();
	double SA(double Frozen_T);
	double SA_new(double Frozen_T, bool flag, size_t Movement);
	double std_dvi();
	double SA_modified(double *rhs);
	Node *find_VDD_spot(Node * nd);
	void form_nbr_pads(Node *rm_pad, vector<Node*>&nbr_pads);
	void update_queue(CircularQueue &q, Node *nd, size_t iter);
	void update_queue_optimize(queue<Node*> &q, Node *nd, size_t iter);
	double update_node_value(int iter, Node *&rm_pad, Node *&nd);
	void update_pad_value(Node *rm_pad, Node *add_pad, 
	vector<Node*>&nodesUpdate_move, int iter_move);
	void update_pad_value_optimize(Node *rm_pad, Node *add_pad, 
		int iter_move, double *rhs);
	void one_move(vector<Node*>&nodesUpdate_move,
	  Node *&rm_pad, Node *&add_pad, 
	  size_t &rm_pad_index, size_t iter_move,
	  Net *&rm_net, Net *&add_net, size_t &index_rm_net);
	void one_move_new(vector<Node*>&nodesUpdate_move,
	  Node *&rm_pad, Node *&add_pad, 
	  size_t &rm_pad_index, size_t iter_move);
	void one_move_modified(vector<Node*>&nodesUpdate_move,
	  double *rhs, Node *&rm_pad, Node *&add_pad,  
	  size_t &rm_pad_index, size_t iter_move);
	double update_cost(vector<Node*> &nodesUpdate_move, 
		int iter_T, double &change_cost_total,
		double *new_voltages);
	double update_cost_modified(vector<Node*> &nodesUpdate_move, double *old_voltages);
	void accept_move(vector<Node*>&nodesUpdate_move, 
		double *new_voltages, size_t rm_index,
		Node *add_pad, Net *&rm_net);
	void accept_move_new(vector<Node*>&nodesUpdate_move, 
		double *new_voltages, size_t rm_index,
		Node *add_pad);
	void reject_move(vector<Node*>&nodesUpdate_move, 
		Node *rm_pad, Node *add_pad, double *new_voltages,
		Net *&rm_net, Net *&add_net, size_t index_rm_net);
	void reject_move_new(vector<Node*>&nodesUpdate_move, 
		Node *rm_pad, Node *add_pad, double *new_voltages);
	void recompute_worst_IRdrop();
	////////////////// end of function for Pad relocation ////
	
	cholmod_common c, *cm;
	size_t peak_mem;
	size_t CK_mem;

	// member functions
	void solve_LU();
	void solve_LU_core();
private:

	bool solve_IT();
	void solve_block_LU();

	bool solve_pcg();
	//bool solve_block_pcg();

public:	
	// initialize things before solve_iteration
	void solve_init();
private:
	// updates nodes value in each iteration
	double solve_iteration();
	void block_init();
	void update_block_geometry();

	// methods of stamping the matrix
	void stamp_by_set(Matrix & A, double* b);
	void stamp_resistor(Matrix & A, Net * net);
	void stamp_current(double* b, Net * net);
	void stamp_VDD(Matrix & A, double* b, Net * net);
	
	void make_A_symmetric(Matrix &A, double *bp);
	void make_A_symmetric_block();

	void stamp_block_matrix();
	void stamp_boundary_matrix();
	void stamp_boundary_net(Net * net);
	void stamp_block_resistor(Net *net, Matrix * A);
	void stamp_block_current(Net * net, Matrix * A);
	void stamp_block_VDD(Net * net, Matrix * A);

	void update_block_rhs(Block & block, int dir);

	//  ******* method for PCG method  ********
	// solve circuit with preconditioned pcg method
	void copy_node_voltages_block(bool from=true);

	// after solving, copy node voltage from replist to nodes
	void get_voltages_from_LU_sol(double* x);
	void get_voltages_from_block_LU_sol();
	void get_vol_mergelist();

	Vec compute_precondition(const Matrix & ML, const Matrix & D, 
			const Matrix & MU, Vec &r);
	void init_precondition(const Matrix &A, Matrix &ML, Matrix &D,
			Matrix &MU);
	// searching all the blocks for global residue
	Vec get_block_res(const Vec& b, const Vec& xk1);
	// searching all the blocks for global zk1
	Vec compute_block_precondition( Vec &r);

	void set_len_per_block();
	void find_block_size ();
	void block_boundary_insert_net(Net * net);
	void find_block_base();

	void partition_circuit();
	double modify_voltage(Block & block, double* x_old);

	void node_voltage_init();
	void solve_one_block(size_t block_id);

	void select_omega();

	void set_type(CIRCUIT_TYPE type){circuit_type = type;};

	void get_samples();

	bool check_diverge() const;

	void merge_along_dir(Node *, DIRECTION dir);
	Node * merge_along_dir_one_pass(Node *, DIRECTION dir, bool remove);
	void merge_node(Node * node);

	// ************** member variables *******************
public:
	NodePtrVector nodelist;		// a set of nodes
private:
	NodePtrVector replist;		// a set of representative nodes
	NodePtrVector mergelist;	// nodes for merging
public:
	NetList net_set[NUM_NET_TYPE];// should be the same as size of NET_TYPE
private:
	// defines the net direction in layers
	static vector<LAYER_DIR> layer_dir;
	vector<int> layers;
	
	// mapping from name to Node object pointer
	unordered_map<string, Node*> map_node;

	// mapping from Net pointer to their index in netlist
	unordered_map<Net*, size_t> net_id;
	unordered_map<Net*, size_t> net_id_pad;

	// circuit name
	string name;

	// blocks
	BlockInfo block_info;
	size_t x_min, y_min, x_max, y_max;

	// control variables
	static double EPSILON;
	static double OMEGA;
	static double OVERLAP_RATIO;
	static size_t MAX_BLOCK_NODES;
	static int MODE; // 0 = IT, 1 = LU

	CIRCUIT_TYPE circuit_type;

	NodePtrVector sample;
public:
	double VDD;
};

inline size_t Circuit::get_total_num_layer(){return layer_dir.size();}

// adds a node into nodelist
inline bool Circuit::add_node(Node * node){
	nodelist.push_back(node);
	map_node[node->name] = node;
	return true;
}

// adds a net into netset
inline bool Circuit::add_net(Net * net){
	if( net->type == RESISTOR )
		net_id[net] = net_set[net->type].size();
	net_set[net->type].push_back(net);
	return true;
}

// adds a net into netset
inline bool Circuit::add_net_pad(Net * net){
	if( net->type == PAD )
		net_id_pad[net] = net_set[net->type].size();
	net_set[net->type].push_back(net);
	return true;
}

// fina a node by name
inline bool Circuit::has_node(string name) const{
	if( map_node.find(name) != map_node.end() ) return true;
	return false;
}

// get a node by name
inline Node * Circuit::get_node(string name){
	unordered_map<string, Node*>::const_iterator it = map_node.find(name);
	if( it != map_node.end() ) return it->second;
	else return NULL;
}

inline void Circuit::merge_node(Node * node){
	for(DIRECTION dir = WEST; dir <= NORTH; dir=DIRECTION(dir+1)){
		// test whether this line has been processed
		if( node->end[dir] != node ) continue;

		// probe for one step, if the line is only one step, don't do it.
		Node * next = node->get_nbr_node(dir);
		if( next == NULL || !next->is_mergeable() ) continue;
		merge_along_dir(node, dir);
	}
}

/*
// find a net by name
inline bool Circuit::has_net(string name) const{
	if( map_net.find(name) != map_net.end() ) return true;
	return false;
}


// get a net by name
inline Net * Circuit::get_net(string name){return map_net[name];}
*/

bool compare_node_ptr(const Node *a, const Node *b);
ostream & operator << (ostream & os, const NodePtrVector & nodelist);
ostream & operator << (ostream & os, const NetList & nets);
//ostream & operator << (ostream & os, const vector<Block > & block_info);

int random_gen(int min, int max);


#endif
