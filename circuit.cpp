// ----------------------------------------------------------------//
// Filename : circuit.cpp
// Author : Zigang Xiao <zxiao2@illinois.edu>
//          Ting Yu <tingyu1@illinois.edu>
//
// implementation file of circuit.h
// ----------------------------------------------------------------//
// - Zigang Xiao - Sun Jan 30 18:35:56 CST 2011
//   * Add UMFPACK support
// - Zigang Xiao - Tue Jan 25 17:19:21 CST 2011
//   * added framework of PCG
// - Zigang Xiao - Tue Jan 18 21:46:13 CST 2011
//   * added solve() and related function
// - Zigang Xiao - Sun Jan 16 16:04:03 CST 2011
//   * added this log

#include <fstream>
#include <string>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <ctime>
#include <utility>
#include <cassert>
#include <vector>
#include "cholmod.h"
#include "umfpack.h"
#include "circuit.h"
#include "util.h"
#include "algebra.h"
#include "node.h"
using namespace std;

double Circuit::EPSILON = 1e-5;
size_t Circuit::MAX_BLOCK_NODES =5500;
double Circuit::OMEGA = 1.2;
double Circuit::OVERLAP_RATIO = 0;//0.2;
int    Circuit::MODE = 0;
const int MAX_ITERATION = 1;//1000;//1000000;//1000;
const int SAMPLE_INTERVAL = 5;
const size_t SAMPLE_NUM_NODE = 10;
const double MERGE_RATIO = 0.3;

//////////////////////////////////////////////////////////////////////////
// Constructor and utility functions goes here

vector<LAYER_DIR> Circuit::layer_dir(MAX_LAYER);

// constructor of Circuit class, name is optional
Circuit::Circuit(string _name):name(_name),
	x_min(INFTY),y_min(INFTY),x_max(0),y_max(0),
	circuit_type(UNKNOWN), VDD(0.0){
	// add ground node
	Node * gnd = new Node(string("0"), Point(-1,-1,-1));
	gnd->rep = gnd;
	this->add_node(gnd);

	for(int i=0;i<MAX_LAYER;i++)
		layer_dir[i]=NA;
	VDD_set.resize(0);
	VDD_candi_set.resize(0);
	CriticalNodes.clear();
	peak_mem = 0;
	CK_mem = 0;
}

// Trick: do not release memory to increase runtime
Circuit::~Circuit(){
	for(size_t i=0;i<nodelist.size();i++) delete nodelist[i];
	for(int type=0;type<NUM_NET_TYPE;type++){
		NetList & ns = net_set[type];
		NetList::iterator it;
		for(it=ns.begin();it!=ns.end();++it)
			delete *it;
	}
	VDD_set.clear();
	VDD_candi_set.clear();
	CriticalNodes.clear();
}

void Circuit::check_sys() const{
	clog<<"**** CHECKING SYSTEM ENVIRONMENT ****"<<endl;
	clog<<"* int size     = "<< sizeof(int)<<endl;
	clog<<"* long size    = "<< sizeof(long)<<endl;
	clog<<"* size_t size  = "<< sizeof(size_t)<<endl;
	clog<<"* UF_long size = "<< sizeof(UF_long)<<endl;
	clog<<"* Max nodelist = "<<(size_t)nodelist.max_size()<<endl;
	clog<<"****            END              ****"<<endl<<endl;
}

// functor to be used in STL sort function
// order: y > x > z > flag 
// input: two node a, b
// return true if a < b, false o/w
// note that ground note are put to last
bool compare_node_ptr(const Node * a, const Node * b){
	if( a->is_ground() ) return false;
	if (b->is_ground() ) return true;

	if( a->pt.y == b->pt.y ){
		if( a->pt.x == b->pt.x ){
			if( a->pt.z == b->pt.z ){
				return a->isX();
			}
			else{
				return (a->pt.z > b->pt.z);// top down
			}
		}
		else
			return ( a->pt.x < b->pt.x );
	}
	else
		return (a->pt.y < b->pt.y);
}

// sort the nodes according to their coordinate 
void Circuit::sort_nodes(){
	sort(nodelist.begin(), nodelist.end(), compare_node_ptr);
	// update node id mapping, 
	// NOTE: ground node will be the last
}

string Circuit::get_name() const{return this->name;}

ostream & operator << (ostream & os, const NodePtrVector & nodelist){
	for(size_t i=0;i<nodelist.size();i++)
		os<<*nodelist[i]<<endl;
	return os;
}

ostream & operator << (ostream & os, const NetList & nets){
	NetList::const_iterator it;
	for(it=nets.begin();it!=nets.end();++it)
		if( (*it) != NULL ) os<<**it<<endl;
	return os;
}

ostream & operator << (ostream & os, const Circuit & ckt){
	os<<"Circuit ["<<ckt.name<<"] info:"<<endl;

	os<<"==== Nodes ===="<<endl;
	os<<ckt.nodelist;

	os<<"==== Reps  ===="<<endl;
	os<<ckt.replist;

	os<<"==== Nets  ===="<<endl;
	os<<ckt.net_set[RESISTOR];

	return os;
}

void Circuit::print(){
	// uncomment this if want to output to a file
	//freopen("output.txt","w",stdout);

	// don't output ground node
	for(size_t i=0;i<nodelist.size()-1;i++){
		printf("%s  %.5e\n", nodelist[i]->name.c_str(), 
				nodelist[i]->value);
	}
}

void Circuit::print_power(){
	// uncomment this if want to output to a file
	//freopen("output.txt","w",stdout);

	// don't output ground node
	for(size_t i=0;i<nodelist.size()-1;i++){
		printf("%s  %.5e\n", nodelist[i]->name.c_str(), 
				nodelist[i]->power);
	}
}

void Circuit::print_matlab(){
	// uncomment this if want to output to a file
	//freopen("output.txt","w",stdout);

	// don't output ground node
	for(size_t i=0;i<nodelist.size()-1;i++){
		printf("%ld %ld  %.5e\n", nodelist[i]->pt.y+1, nodelist[i]->pt.x+1, 
				1-nodelist[i]->value);
	}
}

// find the max IRdrop and min IRdrop numbers, as well as steps
void Circuit::locate_maxIRdrop(){
	max_IRdrop=0;
	double IRdrop = 0;
	for(size_t i=0;i<nodelist.size()-1;i++){
		IRdrop = VDD-nodelist[i]->value;
		if(IRdrop>max_IRdrop)
			max_IRdrop = IRdrop;
	}
}

void Circuit::locate_thIRdrop(){
	th_IRdrop = max_IRdrop*0.8;
}

void Circuit::build_criticalNodes(){
	size_t N = nodelist.size()-1;
	double IRdrop;
	for(size_t i=0;i<N;i++){
		IRdrop = VDD-nodelist[i]->value;
		if(IRdrop > th_IRdrop){
			CriticalNodes.insert(nodelist[i]);
			nodelist[i]->critical = true;
		}
	}
}

// begin the iteration of relaxation
// this is only for 2D case, so there is no nodelist->rep
// for 3D, need to add nodelist->rep
double Circuit::SACost(){
	double cost=0;
	double IRdrop;
	for(size_t i=0;i<nodelist.size();i++){
		IRdrop = VDD-nodelist[i]->value;	
		cost += penalty(IRdrop, max_IRdrop);
	}
	return cost;
}

// simulated annealing
double Circuit::SA(double *rhs){	
	//total cost change of all movement at beginning
	double change_cost_total=0; 
	double P = 0.5; // initial probability
	double T_drop = 0.85; // T=T*T_drop
	// probability to do the movement when cost_change>
	// 0. exp(-cost_change/ Temperature);
	
	// copy the node voltages
	double *new_voltages; 
	new_voltages = new double [nodelist.size()-1];
	for(size_t i=0;i<nodelist.size()-1;i++){
		new_voltages[i] = nodelist[i]->value;
	}
	vector<Node *> nodesUpdate_move;
	Node *rm_pad, *add_pad;
	size_t rm_index;
	rm_pad = NULL;
	add_pad = NULL;	

	double T = 100; // a initial guess
	double Frozen_T=0.01;
	size_t REJECT_LIMIT=30;
	size_t Movement = 30;
	size_t Move_num_rejected=0;
	size_t iter_move=0; size_t iter_T = 0;
	while (T>Frozen_T && Move_num_rejected<REJECT_LIMIT){
		Move_num_rejected = 0;  
		//iter_T;
		for (iter_move=0; iter_move<Movement; iter_move++){
			// compute one movement
			one_move(nodesUpdate_move, rhs, rm_pad, 
				add_pad, rm_index, iter_move);
			double change_cost = update_cost(
				nodesUpdate_move, iter_T, 
				change_cost_total, new_voltages);
			if(change_cost<0)
				accept_move(nodesUpdate_move, 
				  new_voltages, rm_index, add_pad);
			else{
				double prob = exp(-change_cost/T);
				if(acceptProb(prob)==true)
					accept_move(nodesUpdate_move,
					  new_voltages, rm_index, 
					  add_pad);
				else{
					reject_move(nodesUpdate_move,
					  rm_pad, add_pad,
					  new_voltages);
					Move_num_rejected++;
				}
			}
			// recompute worst voltage drop
			recompute_worst_IRdrop(new_voltages);	
		}
		if(iter_T ==0){//calculate the start temperature
			if(change_cost_total >= 0)
				T = -(change_cost_total/Movement)
					/log(P);
			else
				T = (change_cost_total/Movement)
					/log(P);
			//printf("the initial temperature is %f \n", T);
		}//
		//printf("the temperature and probablity of accept is %f, %f \n", T, prob);


		T *= T_drop;
		iter_T++;
	}
	//printf("the total # if iterations: %d \n", iter);
	//printf("the average change of cost is %f \n",  change_cost_total/Movement);
	//printf("the final temperature is %f \n", Temperature);
	//printf("the total temperature changes is %d\n", temp_num);
	//printf("the # of movement at each temperature:%d \n",Movement);

	locate_maxIRdrop();
	locate_thIRdrop();
	double final_cost = SACost();

	nodesUpdate_move.clear();
	delete [] new_voltages;
	return final_cost;
}	
///////////////////////////////////////////////////////////////////////////////
// Computation Functions

// initialization before solving the circuit
// 1. sort the nodes
// 2. set node representatives
// 3. find node in which block, update count
// 4. get representative lists
void Circuit::solve_init(){
	sort_nodes();

	size_t size = nodelist.size() - 1;
	Node * p = NULL;
	for(size_t i=0, nr=0;i<size;i++){
		p=nodelist[i];

		// test if it can be merged
		if( p->is_mergeable() ){
			mergelist.push_back(p);
			continue;
		}

		Net * net = p->nbr[TOP];
		merge_node(p);

		// find the VDD value
		if( p->isX() ) VDD = p->get_value();

		// test short circuit
		if( !p->isX() && // X must be representative 
		    net != NULL &&
		    fzero(net->value) ){
			// TODO: ensure ab[1] is not p itself
			assert( net->ab[1] != p );
			p->rep = net->ab[1]->rep;
		} // else the representative is itself

		// push the representatives into list
		if( p->rep == p ) {
			replist.push_back(p);
			p->rid = nr++;
		}
	}// end of for i

	size_t n_merge = mergelist.size();
	size_t n_nodes = nodelist.size();
	size_t n_reps  = replist.size();
	double ratio = n_merge / (double) (n_merge + n_reps);
	clog<<"mergeable  "<<n_merge<<endl;
	clog<<"replist    "<<n_reps <<endl;
	clog<<"nodelist   "<<n_nodes<<endl;
	clog<<"ratio =    "<<ratio  <<endl;

	net_id.clear();
	net_id_pad.clear();
}

// partition the circuit to X_BLOCKS * Y_BLOCKS blocks
// according to the node size. Assuming they are distributed
// uniformly at random
void Circuit::partition_circuit(){
	size_t num_nodes = replist.size();
	size_t num_blocks =  num_nodes / MAX_BLOCK_NODES;
	if( num_nodes % MAX_BLOCK_NODES > 0 ) ++num_blocks;
	size_t len_x = x_max-x_min;
	size_t len_y = y_max-y_min;
	size_t X_BLOCKS, Y_BLOCKS;

	// Extreme case: only one node in x/y axis
	if( num_blocks == 1 ){
		X_BLOCKS = Y_BLOCKS = 1;
	}
	else if( len_x == 0 ){
		X_BLOCKS = 1;
		Y_BLOCKS = num_blocks;
	}
	else if (len_y == 0 ){
		Y_BLOCKS = 1;
		X_BLOCKS = num_blocks;
	}
	else{// compute length of x and y according to their ratio
		double ratio = double(len_x) / double(len_y);
		X_BLOCKS = sqrt(num_blocks/ratio);
		Y_BLOCKS = X_BLOCKS*ratio;
		if( X_BLOCKS * Y_BLOCKS < num_blocks ) Y_BLOCKS+=1; 
		num_blocks = X_BLOCKS * Y_BLOCKS;
	}
	X_BLOCKS =2;
	Y_BLOCKS =1;
	clog<<"num_blocks: "<<X_BLOCKS<<" / "<<Y_BLOCKS <<endl;
	block_info.X_BLOCKS = X_BLOCKS;
	block_info.Y_BLOCKS = Y_BLOCKS;
	block_info.resize(X_BLOCKS * Y_BLOCKS);
}

// build up block info
// 1. Find block divide point
// 2. Set each node in replist into block
// 3. Compute block size
// 4. Insert boundary netlist into map
void Circuit::block_init(){
	block_info.set_len_per_block(x_min, x_max, y_min, y_max, OVERLAP_RATIO);
	block_info.update_block_geometry();
	find_block_size();
	copy_node_voltages_block();

	//stamp_boundary_matrix();
	stamp_block_matrix();
}

void Circuit::block_boundary_insert_net(Net * net){
	Node * nd[] = {net->ab[0]->rep, net->ab[1]->rep};
	if(nd[0]->is_ground() || nd[1]->is_ground()) return;

	vector<size_t> * ls[] = {&nd[0]->blocklist, &nd[1]->blocklist};

	// idea: for each block where node_a is in
	// find whether node_b is also in the block
	// if yes, then this net is not a boundary net
	// for this block
	vector<size_t>::const_iterator it;
	for(size_t j=0;j<2;j++){
		vector<size_t> *p = ls[j], *q = ls[1-j];
		// for each block k in ls[j]
		for(size_t k=0;k<p->size();k++){
			// find whether block_id in another list
			size_t block_id = (*p)[k];
			it = find( (*q).begin(), (*q).end(), block_id);

			// block_id not in another list, net is boundary
			if( it == (*q).end() ){
				Block & blk = block_info[block_id];
				blk.boundary_netlist.push_back(net);
			}
		}// end of for k
	}// end of for j
}

// stamp the nets by sets, block version
// *NOTE* at the same time insert the net into boundary netlist
void Circuit::stamp_block_matrix(){
	size_t num_blocks = block_info.size();
	Matrix A[num_blocks];
	for(int type=0;type<NUM_NET_TYPE;type++){
		NetList & ns = net_set[type];
		NetList::iterator it;
		switch(type){
		case RESISTOR:
			for(it=ns.begin();it!=ns.end();++it){
				Net * net = *it;
				if( net == NULL ) continue;
				if(net->ab[0]->is_ground() || 
				   net->ab[1]->is_ground()) 
					continue;
				assert( fzero(net->value) == false );
				stamp_block_resistor(*it, A);
			}
			break;
		case CURRENT:
			for(it=ns.begin();it!=ns.end();++it)
				stamp_block_current((*it), A);
			break;
		case VOLTAGE:
			for(it=ns.begin();it!=ns.end();++it){
				if( fzero((*it)->value)  && 
				    !(*it)->ab[0]->is_ground() &&
				    !(*it)->ab[1]->is_ground() )
					continue; // it's a 0v via
				stamp_block_VDD((*it), A);
			}
			break;
		default:
			report_exit("Unknwon net type\n");
			break;
		}
	}
	make_A_symmetric_block();
	clock_t t1, t2;
	t1 = clock();
	// after stamping, convert A to column compressed form
	for(size_t i=0;i<num_blocks;i++){
		if(block_info[i].count>0){
			A[i].set_row(block_info[i].count);		
			block_info[i].CK_decomp(A[i], cm, peak_mem, CK_mem);
		}
	}
	t2 = clock();
	clog<<"decomp time for CK is: "<<1.0*(t2-t1) / CLOCKS_PER_SEC<<endl;
	//clog<<"peak memory for cholmod: "<<peak_mem / 1e9<<" e+06"<<endl;
	//clog<<"CK_mem is: "<<CK_mem / 1e6 <<" G"<<endl;
}

// 1. mark rep nodes into corresponding blocks
// 2. find block size of replist nodes
// 3. allocate rhs size of each block
// 4. find local block index for each node
void Circuit::find_block_size(){
	// start assign replist node into blocks
	size_t n = replist.size();
	Node * p = NULL;
	
	// find which blocks the node belongs to,
	// and its index in the block
	for(size_t i=0;i<n;i++){
		p=replist[i];
		set_blocklist(p);
	}
	// for each block, allocate resource
	for(size_t i=0;i<block_info.size();i++){
		Block & block = block_info[i];
		block.allocate_resource(cm);
		//for(int k=0;k<4;k++) block.A_nbr[k].set_row(block.count);
	}
}

void Circuit::solve(){
	// getting node voltages
	//if( MODE == 0 )
		//solve_IT();
	//else
		solve_LU();

	// compute the power consumption
	// compute_power();
	locate_maxIRdrop();
	locate_thIRdrop();
	// build initial critical nodes list
	build_criticalNodes();
	double cost = 0;
	cost = SACost();
	clog<<"initial cost before SA:" <<cost<<endl;
	double *b;
	b = new double [nodelist.size()-1];
	for(size_t i=0;i<nodelist.size();i++)
		b[i]=0;
	stamp_rhs_SA(b);
	cost = SA(b);
	delete [] b;
	clog<<"final cost after SA:" <<cost<<endl;
}

// restamp current into rhs for SA computation
void Circuit::stamp_rhs_SA(double* b){
	int type=CURRENT;
	NetList & ns = net_set[type];
	NetList::iterator it;

	for(it=ns.begin();it!=ns.end();++it)
		stamp_current(b, (*it));
}

void Circuit::compute_power(){
	NetList::iterator it;
	NetList &ns = net_set[CURRENT];
	for(it = ns.begin(); it!=ns.end();it++){
		Node *nk = (*it)->ab[0]->rep;
		Node *nl = (*it)->ab[1]->rep;
		if(!nk->is_ground()){
			nk->power = nk->value *(*it)->value;
		}
		if(!nl->is_ground()){
			nl->power = nl->value *(*it)->value;
		}
	}
}

// solve Circuit
bool Circuit::solve_IT(){
	// did not find any `X' node
	if( circuit_type == UNKNOWN )
		circuit_type = C4;
	solve_init();

	/*if( replist.size() <= 2*MAX_BLOCK_NODES ){
		clog<<"Replist is small, use direct LU instead."<<endl;
		solve_LU_core();
		return true;
	}
	*/
	select_omega();
	partition_circuit();

	// Only one block, use direct LU instead
	if( block_info.size() == 1 ){
		clog<<"Block size = 1, use direct LU instead."<<endl;
		solve_LU_core();
		return true;
	}

	// cm declared in circuit class
	//cholmod_common c, *cm;
	cm = &c;
	cholmod_start(cm);
	cm->print = 5;
	cm->final_ll = true;
	//cholmod_print_common("first_cm",cm);
	block_init();
	//cholmod_print_common("stamp_cm",cm);

	clog<<"e="<<EPSILON
	    <<"\to="<<OMEGA
	    <<"\tr="<<OVERLAP_RATIO
	    <<"\tb="<<MAX_BLOCK_NODES
	    <<"\tmode="<<MODE<<endl;

	int iter = 0;	
	double diff=0;
	bool successful = false;
	
	clock_t t1, t2;
	t1 = clock();
	while( iter < MAX_ITERATION ){
		diff = solve_iteration();
		iter++;
		clog<<"iter, diff: "<<iter<<" "<<diff<<endl;
		if( diff < EPSILON ){
			successful = true;
			break;
		}
	}
	t2 = clock();
	clog<<"solving iteration use: "<<1.0*(t2-t1)/CLOCKS_PER_SEC<<endl;
	clog<<"# iter: "<<iter<<endl;
	get_voltages_from_block_LU_sol();
	get_vol_mergelist();
	//clog<<"before free. "<<endl;
	for(size_t i=0;i<block_info.size();i++){
		if(block_info[i].count > 0)
			block_info[i].free_block_cholmod(cm);
	}
	//clog<<"after free. "<<endl;
	
	cholmod_finish(cm);
	return successful;
}

// TODO: add comment
void Circuit::node_voltage_init(){
	for(size_t i=0;i<block_info.size();i++){
		Block & block = block_info[i];
		for(size_t j=0;j<block.count;j++){
			block.xp[i] = VDD;
			//block.x[i] = VDD;
			block.nodes[i]->value = VDD;
		}
	}
}

// One iteration during solving the circuit, for any block B:
// 1. update the righthand-side of the matrix of B
// 2. solve the matrix
// 3. update node voltages
// 4. track the maximum error of solution
double Circuit::solve_iteration(){	
	double diff = .0, max_diff = .0;
	for(size_t i=0;i<block_info.size();i++){
		Block &block = block_info[i];
		if( block.count == 0 ) continue;

		block.update_rhs();
		// backup the old voltage value
		
		double *x_old;
		x_old = new double [block.count];
		for(size_t k=0; k<block.count;k++){
			x_old[k] = block.xp[k];
		}
		//cout<<"Matrix A for block: "<<block.bid<<endl;
		block.solve_CK(cm);
		block.xp = static_cast<double *>(block.x_ck->x); 
		//cout<<"block_index: "<<block.bid<<endl;
		//for(size_t j=0;j<block_info[i].count;j++)
			//cout<<j<<" "<<block_info[i].bp[j]<<" "<<block_info[i].xp[j]<<endl;

		// modify node voltage with OMEGA and old voltage value
		diff = modify_voltage(block, x_old);
		delete [] x_old;

		//diff = distance_inf( block.x, x_old );
		if( max_diff < diff ) max_diff = diff;
	}
	return max_diff;
}

double Circuit::modify_voltage(Block & block, double * x_old){
	double max_diff = 0.0;
	OMEGA = 1.0;
	for(size_t i=0;i<block.count;i++){
		block.xp[i] = (1-OMEGA) * x_old[i] + OMEGA * block.xp[i];
		block.nodes[i]->value = block.xp[i];
		double diff = fabs(x_old[i] - block.xp[i]);
		if( diff > max_diff ) max_diff = diff;
	}
	return max_diff;
}

// stamp the matrix and solve
void Circuit::solve_LU_core(){
	size_t n = replist.size();	// replist dosn't contain ground node
	if( n == 0 ) return;		// No node
	//Vec b(n), x(n);
	cholmod_common c, *cm;
	cholmod_dense *b, *x;
	cm = &c;
	cholmod_start(cm);
	cm->print = 5;
	b = cholmod_zeros(n, 1, CHOLMOD_REAL, cm);
	x = cholmod_zeros(n, 1, CHOLMOD_REAL, cm);
	double *bp, *xp;
	bp = static_cast<double *> (b->x);
	Matrix A;
	stamp_by_set(A, bp);
	make_A_symmetric(A, bp);
	//A.merge();

	A.set_row(replist.size());
	Algebra::solve_CK(A, x, b, cm, peak_mem, CK_mem);
	
	xp = static_cast<double *> (x->x);
	// Vec b contains result, copy it back to nodelist
	get_voltages_from_LU_sol(xp);
	get_vol_mergelist();
	cholmod_free_dense(&x, cm);
	cholmod_free_dense(&b, cm);
	cholmod_finish(&c);
}

// solve the node voltages using direct LU
void Circuit::solve_LU(){
	solve_init();
	// build up two VDD pad sets
	pad_set_init();	
	solve_LU_core();
}

// given vector x that obtained from LU, set the value to the corresponding
// node in nodelist
void Circuit::get_voltages_from_LU_sol(double * x){
	for(size_t i=0;i<nodelist.size()-1;i++){
		Node * node = nodelist[i];
		size_t id = node->rep->rid;	// get rep's id in Vec
		double v = x[id];		// get its rep's value
		node->value = v;
	}
}

// compute value of mergelist nodes
void Circuit::get_vol_mergelist(){
	DIRECTION p, q;
	for(size_t i=0;i<mergelist.size();i++){
		Node * node = mergelist[i];
		// check direction
		if( node->nbr[WEST] != NULL ){
			p = WEST;
			q = EAST;
		}
		else{
			p = SOUTH;
			q = NORTH;
		}
		// assign vol value to node
		// left end node and its value
		double r1 = node->eqvr[p];
		double v1 = node->end[p]->value;
		//clog<<" left node: "<<r1<<" / "<<v1<<endl;
		//clog<<" left end "<<node->end[p]->name<<endl;
		// right end node and its value
		double r2 = node->eqvr[q];
		double v2 = node->end[q]->value;
		//clog<<"right node: "<<r2<<" / "<<v2<<endl;
		//clog<<"right end "<<node->end[q]->name<<endl;
		// value for node
		if(v1 > v2){
			node->value = v2 + (v1 - v2) * r2 / (r1 + r2);
		}
		else{
			node->value = v1 + (v2 - v1)  * r1 / (r1 + r2);
		}
		//clog<<" node "<<*node<<endl;
	}
}

// copy solution of block into circuit
void Circuit::get_voltages_from_block_LU_sol(){

	size_t block_id;
	for(size_t i=0;i<nodelist.size()-1;i++){
		Node * node = nodelist[i];
		if( node->is_mergeable() ) continue;
		block_id = node->rep->blocklist[0];
		Block &block = block_info[block_id];
		size_t id = node->rep->id_in_block[0];
		//Vec &p = block.x;
		//double v = p[id];		// get its rep's value
		double v = block.xp[id];
		node->value = v;
	}
}

// solve the circuit using PCG method
/*bool Circuit::solve_pcg(){
	solve_init();
	// remember to exclude ground node
	size_t n = nodelist.size()-1;

	// xk is for k-th iteration, and xk1 is for (k+1)-th iteration
	Vec b(n);
	Vec pk(n), pk1(n);
	Vec rk(n), rk1(n);
	Vec zk(n), zk1(n);
	Vec xk(n), xk1(n);
	double alpha, beta;
	bool successful = false;
	Matrix A, ML, MU, D;

	this->stamp_by_set(A, b);
	
	// initialize first iteration	
	init_precondition(A, ML, D, MU);
	copy_node_voltages(xk1, true);
	rk1 = b - A * xk1; 
	zk1 = compute_precondition(ML, D, MU, rk1);// solve M*z0=r0
	pk1 = zk1;
	int k=0; 
	while(k++ < MAX_ITERATION){
		// new values becomes old values
		zk = zk1;
		pk = pk1;
		xk = xk1;
		double diff=1.0; 

		alpha = (rk * zk) / (pk * A * pk);
		xk1 = xk + alpha * pk;
		rk1 = rk - alpha * (A * pk);
		diff = find_diff(rk1);
		if( diff < EPSILON ) {
			successful = true;
			break;
		}
		zk1 = compute_precondition(ML, D, MU, rk1);
		beta = (rk1 * zk1) / (rk * zk);
		pk1 = zk1 + beta * pk; 
	}
	if( successful ) copy_node_voltages(xk1, false);
	return successful;
}

// M1= (D/OMEGA+L); M2= (D/OMEGA+L'); D:diagonal
void Circuit::init_precondition(const Matrix &A, 
		Matrix &ML, Matrix & D, Matrix &MU){
	A.diagonal_split(ML, D, MU);
	ML *= 1/OMEGA;
	ML += D;

	MU *= 1/OMEGA;
	MU += D;
}

// solve M*z=r; 
// compute a preconditioner for matrix A = ML + D + MU
Vec Circuit::compute_precondition(const Matrix & ML, const Matrix & D, 
			const Matrix & MU, Vec &r){
	size_t n = nodelist.size()-1;
	Vec z1(n), z2(n), z(n);
		
	// solve (D/OMEGA+L)*z=r
	// D is the diagonal of M, L is the lower triangular
	//z1 = Algebra::solve(ML, r);	
	//Algebra::solve(ML, r, z1);	

	// solve (OMEGA/(2-OMEGA)*D^-1)*z=r: D is the diagonal of M,
	z2 = (2-OMEGA)/OMEGA * (D * z1);
	
	// solve (D/OMEGA+L')*z=r: D is the diagonal of M, L' is the upper 
	//z  = Algebra::solve(MU, r);
	//Algebra::solve(MU, r, z);
	return z;
}
*/

// 1. copy node voltages from the circuit to a Vec
//    from = true then copy circuit to x
//    else copy from x to circuit
// 2. map block voltage into global_index
void Circuit::copy_node_voltages_block(bool from){
	size_t id;
	if( from == true ){
		for(size_t i=0;i<replist.size();i++){
			Node *node = replist[i];
			const vector<size_t> &block_id = node->get_block_id();
			for(size_t j=0;j<block_id.size();j++){
				Block &block = block_info[block_id[j]];
				id = node->id_in_block[j];
				block.xp[id] = replist[i]->value;
				//block.x[id] = replist[i]->value;
				block.nodes[id] = replist[i];
			}
		}
	}
	else{
		for(size_t i=0;i<nodelist.size()-1;i++){
			Node * node = nodelist[i];
			const vector<size_t> &block_id = 
				node->rep->get_block_id();
			for(size_t j=0;j<block_id.size();j++){
				Block &block = block_info[block_id[j]];
				id = node->rep->id_in_block[j];
				node->value = block.xp[id];
				//Vec &p = block.x;
				//node->value = p[id];
			}
		}
	}
}

// stamp the net in each set, 
// *NOTE* at the same time insert the net into boundary netlist
void Circuit::stamp_by_set(Matrix & A, double* b){
	for(int type=0;type<NUM_NET_TYPE;type++){
		NetList & ns = net_set[type];
		NetList::iterator it;
		switch(type){
		case RESISTOR:
			for(it=ns.begin();it!=ns.end();++it){
				if( (*it) == NULL ) continue;
				assert( fzero((*it)->value) == false );
				stamp_resistor(A, (*it));
				//block_boundary_insert_net(ns[i]);
			}
			break;
		case CURRENT:
			for(it=ns.begin();it!=ns.end();++it)
				stamp_current(b, (*it));
			break;
		case VOLTAGE:
			for(it=ns.begin();it!=ns.end();++it){
				if( fzero((*it)->value)  && 
				    !(*it)->ab[0]->is_ground() &&
				    !(*it)->ab[1]->is_ground() )
					continue; // it's a 0v via
				stamp_VDD(A, b, (*it));
			}
			break;
		case PAD:
			break;
		default:
			report_exit("Unknwon net type\n");
			break;
		}
	}
}

void Circuit::make_A_symmetric(Matrix &A, double *b){
	int type = RESISTOR;
	NetList & ns = net_set[type];
	NetList::iterator it;
	Node *p, *q;

	for(it=ns.begin();it!=ns.end();it++){
		if( (*it) == NULL ) continue;
			assert( fzero((*it)->value) == false );
		// node a points to X node
		if((*it)->ab[0]->isX()){
			p = (*it)->ab[0]; q = (*it)->ab[1];
		}
		else if((*it)->ab[1]->isX()){
			p = (*it)->ab[1]; q = (*it)->ab[0];
		}
		else continue;
		size_t id = q->rep->rid;
		double G = 1.0 / (*it)->value;
		b[id] += p->value * G;
	}
}

void Circuit::make_A_symmetric_block(){
	int type = RESISTOR;
	NetList & ns = net_set[type];
	NetList::iterator it;
	vector<size_t> *p, *q;
	Node *nk, *nl;

	for(it=ns.begin();it!=ns.end();it++){
		if( (*it) == NULL ) continue;
			assert( fzero((*it)->value) == false );
	
		Net *net = *it;
		Node * nd[] = {net->ab[0]->rep, net->ab[1]->rep};
		if(!nd[0]->isX() && !nd[1]->isX()) continue;
		double G;
		G = 1./net->value;
		vector<size_t> * ls[] = {&nd[0]->blocklist, &nd[1]->blocklist};
		vector<size_t >::const_iterator it_1;
		for(size_t j=0;j<2;j++){
			p = ls[j];
			q = ls[1-j];
			nk = nd[j];
		       	nl = nd[1-j];
			for(size_t i=0;i<p->size();i++){
				// find whether block_id in another list
				size_t block_id = (*p)[i];
				it_1 = find( (*q).begin(), (*q).end(), block_id);
				// 2 nodes in the same block
				if(it_1!=(*q).end() && !nk->isX()){
					Block &block = block_info[block_id];
					size_t k1 = nk->id_in_block[i];
					block.bp[k1] += G *(nl->value);
				}
			}
		}
	}
}
void Circuit::stamp_resistor(Matrix & A, Net * net){
	double G;
	Node * nk = net->ab[0]->rep;
	Node * nl = net->ab[1]->rep;
	size_t k = nk->rid;
	size_t l = nl->rid;
	G = 1./net->value;
	if( !nk->isX() ) {
		A.push_back(k,k, G);
		if(!nl->isX() && l < k) // store lower triangular
			A.push_back(k,l,-G);
	}

	if( !nl->isX() ) {
		A.push_back(l,l, G);
		if(!nk->isX() && k < l) // store ower triangular
			A.push_back(l,k,-G);
	}
}

// stamp a current source
void Circuit::stamp_current(double * b, Net * net){
	Node * nk = net->ab[0]->rep;
	Node * nl = net->ab[1]->rep;
	if( !nk->is_ground() && !nk->isX() ) { 
		size_t k = nk->rid;
		b[k] += -net->value;
	}
	if( !nl->is_ground() && !nl->isX() ) {
		size_t l = nl->rid;
		b[l] +=  net->value;
	}
}

// stamp a voltage source
void Circuit::stamp_VDD(Matrix & A, double * b, Net * net){
	// find the non-ground node
	Node * X = net->ab[0];
	if( X->is_ground() ) X = net->ab[1];
	size_t id = X->rep->rid;
	A.push_back(id, id, 1.0);
	Net * south = X->rep->nbr[SOUTH];
	if( south != NULL &&
	    south->type == CURRENT ){
		// this node connects to a VDD and a current
		assert( feqn(1.0, b[id]) ); // the current should be stamped
		b[id] = net->value;	    // modify it
	}
	else
		b[id] += net->value;
}

// =========== stamp block version of matrix =======

void Circuit::stamp_block_resistor(Net * net, Matrix * A){
	Node * nd[] = {net->ab[0]->rep, net->ab[1]->rep};

	double G;	
	G = 1./net->value;

	vector<size_t> * ls[] = {&nd[0]->blocklist, &nd[1]->blocklist};
	vector<size_t>::const_iterator it;
	for(size_t j=0;j<2;j++){
		vector<size_t> *p = ls[j], *q = ls[1-j];
		Node *nk = nd[j], *nl = nd[1-j];
		for(size_t i=0;i<p->size();i++){
			// find whether block_id in another list
			size_t block_id = (*p)[i];
			it = find( (*q).begin(), (*q).end(), block_id);

			// block_id not in another list, net is boundary
			if( it == (*q).end() ){
				Block & blk = block_info[block_id];
				blk.boundary_netlist.push_back(net);
				if( !nk->isX() ) {
					// stamp value into block_ids
					size_t k1 = nk->id_in_block[i];
					Matrix &pk = A[block_id];	
					pk.push_back(k1,k1, G);
				}
			}
			// else 2 nodes belongs to the same block
			// stamp resistor
			else if( !nk->isX() ) {
				size_t k1 = nk->id_in_block[i];
				Matrix &pk = A[block_id];

				size_t j1 = it - (*q).begin();
				size_t l1 = nl->id_in_block[j1];

				pk.push_back(k1,k1, G);
				if(!nl->isX() && l1 < k1) // only store the lower triangular part
					pk.push_back(k1,l1,-G);
			}
		}// end of for k
	}// end of for j	
}

void Circuit::stamp_block_current(Net * net, Matrix * A){
	Node * nk = net->ab[0]->rep;
	Node * nl = net->ab[1]->rep;
	if( !nk->is_ground() && !nk->isX() ) { 
		for(size_t i=0;i<nk->blocklist.size();i++){
			size_t block_idk = nk->blocklist[i];
			Block &block_k = block_info[block_idk];
			//Vec & pk = block_k.b;
			size_t k = nk->id_in_block[i];
			block_k.bp[k] += -net->value;
			//pk[k] += -net->value;
		}
	}
	if( !nl->is_ground() && !nl->isX() ) {
		for(size_t i=0;i<nl->blocklist.size();i++){
			size_t block_idl = nl->blocklist[i];
			Block & block_l = block_info[block_idl];
			//Vec & pl = block_l.b;
			size_t l = nl->id_in_block[i];
			block_l.bp[l] += net->value;
			//pl[l] +=  net->value;
		}
	}
}

void Circuit::stamp_block_VDD(Net * net, Matrix * A){
	// find the non-ground node
	Node * X = net->ab[0];
	if( X->is_ground() ) X = net->ab[1];
	for(size_t i=0;i<X->rep->id_in_block.size();i++){
		size_t block_id = X->rep->blocklist[i];
		Block &block = block_info[block_id];	
		Matrix & p = A[block_id];
		//Vec & q = block.b;
		size_t id =X->rep->id_in_block[i];
		p.push_back(id, id, 1.0);
		Net * south = X->rep->nbr[SOUTH];
		if( south != NULL &&
	   	 south->type == CURRENT ){
		// this node connects to a VDD and a current
		// the current should be stamped
			//assert( feqn(1.0, q[id]) ); 
			assert( feqn(1.0, block.bp[id]) );
			block.bp[id] = net->value;
			//q[id] = net->value;	    // modify it
		}
		else{
			block.bp[id] += net->value;
			//q[id] += net->value;
		}
	}
}

// set the block_id in a node
// according to the following order
// where 0 is the original block the node should be 
// 8 1 2
// 7 0 3
// 6 5 4
//
void Circuit::set_blocklist(Node * nd){
	const double len_per_block_x = block_info.len_per_block_x;
	const double len_per_block_y = block_info.len_per_block_y;
	const double len_ovr_x = block_info.len_ovr_x;
	const double len_ovr_y = block_info.len_ovr_y;
	const size_t X_BLOCKS = block_info.X_BLOCKS;
	const size_t Y_BLOCKS = block_info.Y_BLOCKS;
	const long x = nd->pt.x - x_min;	// point relative location
	const long y = nd->pt.y - y_min;
	size_t bx0 = x / len_per_block_x;
	size_t by0 = y / len_per_block_y;
	long bx, by;		// block index
	double lx, ly, ux, uy;	// block bounding box

	const long dx[]={0, 0, 1, 1,  1,  0, -1, -1, -1};
	const long dy[]={0, 1, 1, 0, -1, -1, -1,  0, 1};

	// test whether the node is in one of the nine blocks
	for(int i=0;i<9;i++){
		bx = bx0 + dx[i];
		by = by0 + dy[i];

		// check the block index is valid
		if( bx < 0 || bx >= (long)X_BLOCKS ) continue;
		if( by < 0 || by >= (long)Y_BLOCKS ) continue;

		// compute block coordinate
		lx = bx * len_per_block_x - len_ovr_x;
		ly = by * len_per_block_y - len_ovr_y;
		ux = (bx+1) * len_per_block_x + len_ovr_x;
		uy = (by+1) * len_per_block_y + len_ovr_y;

		// check if the point is in the block
		if( !(x>=lx && x<=ux && y>=ly && y<=uy) ) continue;	

		size_t id = by * X_BLOCKS + bx;
		assert( id<X_BLOCKS* Y_BLOCKS );

		/*
		vector<size_t >::const_iterator it;
		it = find(nd->blocklist.begin(),nd->blocklist.end(), id);
		if(it!=nd->blocklist.end()){
			printf("id=%ld, i=%d\n", id,i);
			printf("xbase,ybase=%lf %lf\n", lx_base, ly_base);
			printf("lx,ly=%lf %lf\n", lx, ly);
			printf("ux,uy=%lf %lf\n", ux, uy);
			printf("x,y=%ld %ld\n", x, y);
			printf("bx,by=%ld %ld\n", bx, by);
			printf("xc,yc=%ld %ld\n", x_center, y_center);
			continue;
		}
		*/
		nd->blocklist.push_back(id);
		Block & block = block_info[id];
		nd->id_in_block.push_back(block.count++);
	}
}

void Circuit::get_parameters(
		double & epsilon,
		double & omega,
		double & overlap_ratio,
		size_t & max_block_nodes,
		int & mode){
	epsilon		= EPSILON;
	omega		= OMEGA;
	overlap_ratio	= OVERLAP_RATIO; 
	max_block_nodes	= MAX_BLOCK_NODES;
	mode		= MODE;
}

// default values of these parameters are at the begining of this file
void Circuit::set_parameters(
		double epsilon, 
		double omega, 
		double overlap_ratio,
		size_t max_block_nodes,
		int mode){
	EPSILON		= epsilon;
	OMEGA		= omega;
	OVERLAP_RATIO 	= overlap_ratio;
	MAX_BLOCK_NODES	= max_block_nodes;
	MODE		= mode;
}

// choose an appropriate omega for the circuit s.t.
// - node size (use replist)
// - type (c4 or wb)
// - number of layers
void Circuit::select_omega(){
	double omega=OMEGA;
	size_t num_nodes = replist.size();
	size_t num_layers = layers.size();
	if( num_nodes < 0.05e6 )
		omega=1.0;
	else if (num_nodes < 0.2e6 )
		omega = 1.1;
	else if (num_nodes < 0.3e6 )
		omega = 1.2;
	else if (num_nodes < 0.5e6 )
		omega = 1.3;
	else if (num_nodes < 1.2e6 )
		omega = 1.4;
	else
		omega = 1.5;

	if( circuit_type == WB && num_layers >= 8 ) omega += 0.2;

	if( circuit_type == C4 ) omega += 0.1;

	if( name == "GND" && num_nodes < 1.2e6) omega -= 0.1;

	if( omega >= 1.6 ) omega = 1.6;
	if( omega <= 1.0 ) omega = 1.0;

	OMEGA = omega;
}

// Randomly choose a number of sample nodes to monitor
void Circuit::get_samples(){
	size_t num_nodes = replist.size();
	srand(time(NULL));
	while(sample.size()<SAMPLE_NUM_NODE){
		int id = rand() % num_nodes;
		sample.push_back(replist[id]);
	}
}

bool Circuit::check_diverge() const{
	for(size_t i=0;i<SAMPLE_NUM_NODE;i++){
		double x = sample[i]->value;
		if(VDD > 0){
			if( x < 0.0 || x > VDD ) return true;
		}
		else
			if(x<0.0) return true;
	}
	return false;
}

Node * Circuit::merge_along_dir_one_pass(Node * start, DIRECTION dir, bool remove){
	double sum = 0.0;
	DIRECTION ops = get_opposite_dir(dir);
	Node * p = start;

	// traverse along the direction, sum the resistor value and set the node end
	while(1){
		p = p->get_nbr_node(dir);
		p->end[ops] = start;
		Net * net = p->nbr[ops];
		sum += net->value;
		p->eqvr[ops] = sum;
		if( remove ) {
			size_t id = net_id[net];
			net_set[RESISTOR][id] = NULL;
			delete net;
		}
		if( !p->is_mergeable() ) break;
	}

	return p;	// return end point
}

// merge a line along direction
void Circuit::merge_along_dir(Node * node, DIRECTION dir){
	// two pass traversal
	DIRECTION ops = get_opposite_dir(dir);
	node->end[dir] = merge_along_dir_one_pass(node, dir, false);
	Node * other = node->end[dir];
	other->end[ops] = node;
	merge_along_dir_one_pass(other, ops, true);
	//assert( ret == node );

	// add a new net between `node' and its end
	Net * net = new Net(RESISTOR, node->eqvr[dir], node, other);
	node->nbr[dir] = other->nbr[ops] = net;
	this->add_net(net);
}

inline double Circuit::penalty(double v, double vworst){
	double penalty;
	if(v>=0.8*vworst && v<0.98*vworst)
		penalty = v*v*10;
	else if(v>=0.98*vworst)
		penalty = v*v*100;
	else
		penalty = 0;
	return penalty;
}

// build up VDD_set and VDD_candi_set
void Circuit::pad_set_init(){
	for(size_t i=0;i<nodelist.size();i++){
		if(nodelist[i]->is_candi()==true)
			VDD_candi_set.push_back(nodelist[i]);
		if(nodelist[i]->isX()==true)
			VDD_set.push_back(nodelist[i]);
	}	
}

void Circuit::one_move(vector<Node*>&nodesUpdate_move,
	double *rhs, Node *rm_pad, Node *add_pad,  
	size_t &rm_pad_index, size_t iter_move){
	size_t VDD_num = VDD_set.size();
	vector<Node*> nbr_pads;
	size_t nbr_index;
	Node *ref_node;

	// the nodes updated in each movement 
	nodesUpdate_move.resize(0);

	//1. randomly pick a Vdd pad to remove
	rm_pad_index = random_gen(0, VDD_num-1); 
	rm_pad = VDD_set[rm_pad_index];	
	// It is no more X node
	rm_pad->flag= false; 

	// find its neiboring pad candidates	
	form_nbr_pads(rm_pad, nbr_pads);	
	if(nbr_pads.size()==0) return;
	
	// random select a valid nbr pad
	nbr_index = random_gen(0, nbr_pads.size()-1);
	add_pad = nbr_pads[nbr_index];
	add_pad->flag = true; 
	rm_pad->value=0;
	add_pad->value=VDD;	
	ref_node = rm_pad;
	// update pad value and nbr area by iterations
	update_pad_value(rm_pad, add_pad, nodesUpdate_move, 
			iter_move, rhs);
}

void Circuit::form_nbr_pads(Node *rm_pad, vector<Node*>&nbr_pads){
	Node *na, *nb;
	Net *net;
	//find its neighboring pad nodes
	nbr_pads.resize(0);
	for(size_t j=0;j<4;j++){
		if(rm_pad->nbr_pad[j]!=NULL){
			net = rm_pad->nbr_pad[j];
			na = net->ab[0];
			nb = net->ab[1];
			if(rm_pad->name==na->name &&
					nb->isX()==false)
				nbr_pads.push_back(nb);
			else if(rm_pad->name==
					nb->name && na->isX()==false)
				nbr_pads.push_back(na);	
		}
	}
}

void Circuit::update_pad_value(Node *rm_pad, Node *add_pad, 
	vector<Node*>&nodesUpdate_move, int iter_move, double *rhs){
	size_t iter = 0;
	double eps0=1e-4;
	double eps1 = 1e-4;
	size_t LIMIT = 500;
	double V_refer=0;
	double V_refer_old=0;
	Node *ref_node;
	ref_node = rm_pad;
	
	int *timestamps;
	timestamps = new int [nodelist.size()-1];
	for(size_t i=0;i<nodelist.size();i++)
		timestamps[i]=0;

	size_t MAX_QUEUE_SIZE = nodelist.size()-1;
	CircularQueue q(MAX_QUEUE_SIZE);

	while(fabs(V_refer-ref_node->value) > eps0 && iter<LIMIT){
		V_refer_old = V_refer;
		V_refer = ref_node->value;
		q.reset();
		q.insert(rm_pad);
		q.insert(add_pad);
		rm_pad->flag_visited = iter; 
		add_pad->flag_visited = iter;
		double V_improve = 0;
		while(!q.isEmpty() && (V_improve>eps1)){
			Node *nd = q.extractFront();
			if(timestamps[nd->rid]!=iter_move){
				nodesUpdate_move.push_back(nd);
				timestamps[nd->rid] = iter_move;
			}
			// update node value with neighbors
			V_improve = update_node_value(iter, 
					rm_pad, nd, rhs);
			// update_queue with nd's neighbors
			update_queue(q, nd, iter);	
		}
		iter++;
	}
	delete [] timestamps;	
}

// rhs is the current vector
double Circuit::update_node_value(int iter, Node *rm_pad, Node *nd, double *rhs){
	if(nd->isX()==true) return 0;
	double h = 0.06;
	double omega = 2-h;
	printf("the omega and h: %3.10lf %3.10lf\n",omega,h);
	double V_old=0;
	double V_temp = 0;
	double G = 0;
	Net *net;
	Node *nbr, *na, *nb;
	double sum = 0;
	
	V_old = nd->value; 
	// update nd->value
	for(int i=0;i<6;i++){
		net = nd->nbr[i];
		if(net ==NULL) continue;
		G = 1.0/net->value;
		sum += G;
		na = net->ab[0]; nb = net->ab[1];
		if(nd->name == na->name) nbr = nb;
		else	nbr = na;
		V_temp += G*nbr->value;
	}
	V_temp += rhs[nd->rid];
	V_temp /=sum;
	if(iter ==0 && nd->name==rm_pad->name)
		nd->value  = V_temp;
	else
		nd->value = (1-omega)*nd->value + omega*V_temp;
 	double V_improve = fabs(nd->value - V_old);
	return V_improve;
}

void Circuit::update_queue(CircularQueue &q, Node *nd, size_t iter){
	Net * net; Node *nbr;
	Node *na, *nb;
	for(int i=0;i<4;i++){
		net = nd->nbr_pad[i];
		if(net==NULL) continue;
		// find the other node, which is nbr
		na = net->ab[0];
		nb = net->ab[1];
		if(nd->name == na->name)
			nbr = nb;
		else	nbr = na;
		if(nbr->flag_visited != iter){
			q.insert(nbr);
			nbr->flag_visited = iter;
		}
	}
}

// update cost and total_cost
double Circuit::update_cost(vector<Node*> &nodesUpdate_move, 
		int iter_T, double &change_cost_total,
		double *new_voltages){
	double change_cost = 0;
	Node *nd;
	for (size_t i=0; i<nodesUpdate_move.size();i++){
		nd = nodesUpdate_move[i] ;
		change_cost += penalty(VDD-nd->value, max_IRdrop)-
			penalty(VDD-new_voltages[nd->rid], max_IRdrop);
	}
	if (iter_T==0){
		// accept all the movements at the beginning
		change_cost_total += change_cost;
	}
	return change_cost;
}

void Circuit::accept_move(vector<Node*>&nodesUpdate_move, 
		double *new_voltages, size_t rm_index, 
		Node *add_pad){
	Node *nd;
	for (size_t i=0; i<nodesUpdate_move.size();i++){
		nd = nodesUpdate_move[i] ;

		if (new_voltages[nd->rid] != nd->value) {
			if (nd->critical==true){
				// remove the node from critical 
				// set O(logN)
				CriticalNodes.erase(nd) ; 
				nd->critical = false ; 
			}
			new_voltages[nd->rid] = nd->value; 
			if ((VDD - nd->value) > th_IRdrop) {
				// insert it into the critical set 
				// O(logN)
				CriticalNodes.insert(nd) ; 
				nd->critical = true ;
			}
		}

	}
	//keep the change of Vdd
	VDD_set[rm_index]=add_pad;
	   
        //printf("movement %d accepted because change of cost 
        //%f\n", move, change_cost);
}

void Circuit::reject_move(vector<Node*>&nodesUpdate_move, 
	Node *rm_pad, Node *add_pad, double *new_voltages){
	Node *nd;
	for (size_t i=0; i<nodesUpdate_move.size();i++){
		nd = nodesUpdate_move[i];
		nd->value = new_voltages[nd->rid];
	}
	rm_pad->flag = true; // rm_pad isX again
	add_pad->flag = false; //add_pad is not X again
	//printf("movement %d rejected because change of cost %f\n",
	//move, change_cost); 
}

void Circuit::recompute_worst_IRdrop(double *new_voltages){
	Node *nd;
	if (!CriticalNodes.empty()) {
		nd = *CriticalNodes.begin() ; // the first element in criticalNodes will be the node with smallest volt
		max_IRdrop = new_voltages[nd->rid] ;
	}
	else {
		locate_maxIRdrop();
		locate_thIRdrop();
		build_criticalNodes();	
	} 
	//printf("the worst voltage drop: %f\n", max_IRdrop);
}

bool Circuit::acceptProb(double p){
	return rand()<p*((double)RAND_MAX+1.0);
}

int random_gen(int min, int max){
    return min + int( ((max-min) +1) * (rand() /(RAND_MAX + 1.0))) ; 

}
