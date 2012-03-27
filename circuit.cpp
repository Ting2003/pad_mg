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
	map_node.clear();
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
double Circuit::locate_maxIRdrop(){
	max_IRdrop=0;
	double IRdrop = 0;
	double sum = 0;
	for(size_t i=0;i<nodelist.size()-1;i++){
		//clog<<"i, value: "<<i<<" "<<*nodelist[i]<<endl;
		IRdrop = fabs(VDD-nodelist[i]->value);
		sum += IRdrop;
		//clog<<"i, cur, max: "<<i<<" "<<IRdrop<<" "<<max_IRdrop<<endl;
		if(IRdrop>max_IRdrop)
			max_IRdrop = IRdrop;
	}
	return sum;
	//clog<<"max IRdrop "<<max_IRdrop<<endl;
}

void Circuit::locate_thIRdrop(){
	th_IRdrop = max_IRdrop*0.8;
}

void Circuit::build_criticalNodes(){
	CriticalNodes.clear();
	size_t N = nodelist.size()-1;
	double IRdrop;
	for(size_t i=0;i<N;i++){
		IRdrop = VDD-nodelist[i]->value;
		if(IRdrop > th_IRdrop){
			CriticalNodes.insert(nodelist[i]);
			nodelist[i]->critical = true;
		}
		else nodelist[i]->critical = false;
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
double Circuit::SA(double Frozen_T){	
	//total cost change of all movement at beginning
	double change_cost_total=0; 
	double P = 0.5; // initial probability
	double T_drop = 0.85; // T=T*T_drop
	// probability to do the movement when cost_change>
	// 0. exp(-cost_change/ Temperature);
	
	// copy the node voltages
	size_t N = nodelist.size()-1;
	double *old_voltages; 
	old_voltages = new double [N];
	for(size_t i=0;i<N;i++){
		old_voltages[i] = nodelist[i]->value;
	}
	vector<Node *> nodesUpdate_move;
	Node *rm_pad, *add_pad;
	size_t rm_index;
	rm_pad = NULL;
	add_pad = NULL;	

	double T = 100; // a initial guess
	//double Frozen_T=0.001;
	size_t Movement = 10;
	size_t Move_num_rejected=0;
	size_t iter_move=1; size_t iter_T = 0;
	double prev_maxIRdrop = max_IRdrop;
	double prev_maxIRdrop_T = max_IRdrop;

	double prev_sumIRdrop = locate_maxIRdrop();
	double cur_sumIRdrop = 0;
	
	
	//clog<<"before starting T iteration."<<endl;
	while (T > Frozen_T){
		Move_num_rejected = 0;  
		for (iter_move=1; iter_move<Movement; iter_move++){
			// compute one movement
			one_move(nodesUpdate_move, rm_pad, 
				add_pad, rm_index, iter_move);
			cur_sumIRdrop = locate_maxIRdrop();
			//clog<<"new max: "<<max_IRdrop<<endl;
			double change_cost = max_IRdrop - prev_maxIRdrop;
			double change_cost_total = cur_sumIRdrop - 
				prev_sumIRdrop;
			//double change_cost = update_cost(
				//nodesUpdate_move, iter_T, 
				//change_cost_total, old_voltages);
			if(change_cost<0 || change_cost_total < 0){
				//clog<<"accept. "<<endl;
				if(change_cost <0)
					prev_maxIRdrop = max_IRdrop;
				if(change_cost_total < 0)
					prev_sumIRdrop = cur_sumIRdrop;
				accept_move(nodesUpdate_move, 
				  old_voltages, rm_index, add_pad);
			}
			else{
				double prob = exp(-change_cost/T);
				if(acceptProb(prob)==true){								accept_move(nodesUpdate_move,
					  old_voltages, rm_index, 
					  add_pad);
					  //locate_maxIRdrop();
					  //if(max_IRdrop < prev_maxIRdrop)
					  //prev_maxIRdrop = max_IRdrop;
				}else{
					reject_move(nodesUpdate_move,
					  rm_pad, add_pad,
					  old_voltages);
					Move_num_rejected++;
				}
			}
			// recompute worst voltage drop
			//recompute_worst_IRdrop();
		}
		//clog<<"change_cost_total: "<<change_cost_total<<endl;
		/*if(iter_T ==1){//calculate the start temperature
			if(change_cost_total >= 0)
				T = -(change_cost_total/Movement)
					/log(P);
			else
				T = (change_cost_total/Movement)
					/log(P);
			//printf("the initial temperature is %f \n", T);
		}//*/
		//printf("the temperature and probablity of accept is %f, %f \n", T, prob);

		T *= T_drop;
		//clog<<endl<<"origin max, origin_total: "<<prev_maxIRdrop<<endl;
		//clog<<"iter_T, T, stop_prob: "<<iter_T<<" "<<T<<" "<<Move_num_rejected<<" / "<<Movement<<endl;
		if(1.0*Move_num_rejected / Movement >= 0.99) break;
		iter_T++;
	}

	//locate_maxIRdrop();
	//locate_thIRdrop();
	//double final_cost = SACost();
	//clog<<"final_cost: "<<max_IRdrop<<endl;

	nodesUpdate_move.clear();
	delete [] old_voltages;
	return max_IRdrop;
}	
// simulated annealing
double Circuit::SA_new(double Frozen_T){	
	//total cost change of all movement at beginning
	double change_cost_total=0; 
	double P = 0.5; // initial probability
	double T_drop = 0.85; // T=T*T_drop
	// probability to do the movement when cost_change>
	// 0. exp(-cost_change/ Temperature);
	
	// copy the node voltages
	size_t N = nodelist.size()-1;
	double *old_voltages; 
	old_voltages = new double [N];
	for(size_t i=0;i<N;i++){
		old_voltages[i] = nodelist[i]->value;
	}
	vector<Node *> nodesUpdate_move;
	Node *rm_pad, *add_pad;
	size_t rm_index;
	rm_pad = NULL;
	add_pad = NULL;	

	double T = 100; // a initial guess
	//double Frozen_T=0.001;
	size_t Movement = 10;
	size_t Move_num_rejected=0;
	size_t iter_move=1; size_t iter_T = 0;
	double prev_maxIRdrop = max_IRdrop;
	double prev_maxIRdrop_T = max_IRdrop;

	double prev_sumIRdrop = locate_maxIRdrop();
	double cur_sumIRdrop = 0;
	
	
	//clog<<"before starting T iteration."<<endl;
	while (T > Frozen_T){
		Move_num_rejected = 0;  
		for (iter_move=1; iter_move<Movement; iter_move++){
			// compute one movement
			one_move(nodesUpdate_move, rm_pad, 
				add_pad, rm_index, iter_move);
			cur_sumIRdrop = locate_maxIRdrop();
			//clog<<"new max: "<<max_IRdrop<<endl;
			double change_cost = max_IRdrop - prev_maxIRdrop;
			double change_cost_total = cur_sumIRdrop - 
				prev_sumIRdrop;
			//double change_cost = update_cost(
				//nodesUpdate_move, iter_T, 
				//change_cost_total, old_voltages);
			if(change_cost<0 || change_cost_total < 0){
				//clog<<"accept. "<<max_IRdrop<<" "<<cur_sumIRdrop<<endl;
				if(change_cost <0)
					prev_maxIRdrop = max_IRdrop;
				if(change_cost_total < 0)
					prev_sumIRdrop = cur_sumIRdrop;
				accept_move(nodesUpdate_move, 
				  old_voltages, rm_index, add_pad);
			}
			else{
				//double prob = exp(-change_cost/T);
				//if(acceptProb(prob)==true){								accept_move(nodesUpdate_move,
					 // old_voltages, rm_index, 
					  //add_pad);
					  //locate_maxIRdrop();
					  //if(max_IRdrop < prev_maxIRdrop)
					  //prev_maxIRdrop = max_IRdrop;
				//}else
				{
					reject_move(nodesUpdate_move,
					  rm_pad, add_pad,
					  old_voltages);
					Move_num_rejected++;
				}
			}
			// recompute worst voltage drop
			//recompute_worst_IRdrop();
		}
		//clog<<"change_cost_total: "<<change_cost_total<<endl;
		/*if(iter_T ==1){//calculate the start temperature
			if(change_cost_total >= 0)
				T = -(change_cost_total/Movement)
					/log(P);
			else
				T = (change_cost_total/Movement)
					/log(P);
			//printf("the initial temperature is %f \n", T);
		}//*/
		//printf("the temperature and probablity of accept is %f, %f \n", T, prob);

		T *= T_drop;
		//clog<<endl<<"origin max, origin_total: "<<prev_maxIRdrop<<" "<<prev_sumIRdrop<<endl;
		//clog<<"iter_T, T, stop_prob: "<<iter_T<<" "<<T<<" "<<Move_num_rejected<<" / "<<Movement<<endl;
		if(1.0*Move_num_rejected / Movement >= 0.9) break;
		iter_T++;
	}

	//locate_maxIRdrop();
	//locate_thIRdrop();
	//double final_cost = SACost();
	//clog<<"final_cost: "<<max_IRdrop<<endl;

	nodesUpdate_move.clear();
	delete [] old_voltages;
	return max_IRdrop;
}	


// simulated annealing
double Circuit::SA_modified(double *rhs){	
	// copy the node voltages
	size_t N = nodelist.size()-1;
	double *old_voltages; 
	old_voltages = new double [N];
	for(size_t i=0;i<N;i++){
		old_voltages[i] = nodelist[i]->value;
	}
	vector<Node *> nodesUpdate_move;
	Node *rm_pad, *add_pad;
	size_t rm_index;
	rm_pad = NULL;
	add_pad = NULL;	

	size_t Movement = 10;
	size_t Move_num_rejected=0;
	size_t iter_move=1;
	bool flag = false;
	
	locate_maxIRdrop();
	locate_thIRdrop();
	build_criticalNodes();
	double cost = 0;
	cost = SACost();
	//cout<<"cost before SA_modified: "<<cost<<endl;
	for(size_t i = 0; i<VDD_set.size();i++){
		//clog<<"i, node: "<<i<<" "<<*VDD_set[i]<<endl;
		// pick up one rm_pad
		//if(VDD_set[i]->flag_visited >=1) continue;
		if(~(VDD_set[i]->flag_visited >=1)){
			rm_pad = VDD_set[i];
			rm_pad->flag = false;
			rm_pad->value = 0;
			rm_index = i;
			clog<<endl<<"rm_pad "<<*rm_pad<<endl;
		}
		Move_num_rejected = 0;  
		for (iter_move=1; iter_move<Movement; iter_move++){
			rm_pad = VDD_set[i];
			rm_pad->flag = false;
			rm_pad->value = 0;
			//clog<<"rm_pad "<<*rm_pad<<endl;

			//clog<<"iter_move: "<<iter_move<<endl;
			if(1.0*Move_num_rejected / Movement >= 0.5){
				flag = true;
				break;
			}
			
			// compute one movement
			one_move_modified(nodesUpdate_move, rhs, 
			  rm_pad, add_pad, rm_index, iter_move);
			//clog<<"one move modified. "<<endl;
			double change_cost = update_cost_modified
			  (nodesUpdate_move, old_voltages);
			//cout<<"move, change_cost: "<<iter_move<<" "<<change_cost<<endl;
			if(change_cost<0){

				//cout<<"rm_pad: "<<*rm_pad<<endl;
				VDD_set[rm_index]=add_pad;
				for(size_t k=0;k<N;k++){
					old_voltages[k] = nodelist[k]->value;
				}
				//accept_move(nodesUpdate_move, 
				  //old_voltages, rm_index, add_pad);
				locate_maxIRdrop();
				locate_thIRdrop();
				double final_cost = SACost();
				//cout<<"rm_index, new VDD: "<<rm_index<<" "<<*VDD_set[rm_index]<<endl;
				//cout<<"add_pad, change_cost, cost:  "<<*add_pad<<" "<<change_cost<<" "<<final_cost<<endl;
				//clog<<"add_pad: "<<*add_pad<<endl;
			}
			else{
				//clog<<"reject move: "<<endl;
				reject_move(nodesUpdate_move,
				  rm_pad, add_pad, old_voltages);
					Move_num_rejected++;
			}
		}

		// recompute worst voltage drop
		recompute_worst_IRdrop();
		if(flag == true){
			continue;
		}
	}

	locate_maxIRdrop();
	locate_thIRdrop();
	double final_cost = SACost();
	clog<<"final_cost: "<<final_cost<<endl;

	nodesUpdate_move.clear();
	delete [] old_voltages;
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
	//clog<<"decomp time for CK is: "<<1.0*(t2-t1) / CLOCKS_PER_SEC<<endl;
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

void Circuit::solve_coarse(double Frozen_T){
	solve_LU();
	locate_maxIRdrop();
	clog<<"initial max_IRdrop is: 		 "<<max_IRdrop<<endl;
	
	//RANSAC_init();
	//locate_maxIRdrop();
	//clog<<"max_IRdrop after ransac init:	 "<<max_IRdrop<<endl;

	// optimized method plus SA
	optimize_pad_assign_new();
	locate_maxIRdrop();
	clog<<"max_IRdrop after opti:		 "<<max_IRdrop<<endl;
	
	SA(Frozen_T);
	locate_maxIRdrop();
	clog<<"max_IRdrop after SA  :		 "<<max_IRdrop<<endl;
}

void Circuit::solve(double Frozen_T){
	// getting node voltages
	//if( MODE == 0 )
		//solve_IT();
	//else
		solve_LU();
	locate_maxIRdrop();
	//locate_thIRdrop();
	clog<<"initial max_IRdrop is: 	"<<max_IRdrop<<endl;
	// use ransac method to get a better init pad assignment
	/*RANSAC_init();
	locate_maxIRdrop();
	clog<<"max_IRdrop after ransac init:	 "<<max_IRdrop<<endl;

	//Mean_shift_move();
	// optimized method plus SA
	optimize_pad_assign_new();
	locate_maxIRdrop();
	clog<<"max_IRdrop after opti:	"<<max_IRdrop<<endl;
	*/
	SA(Frozen_T);
	locate_maxIRdrop();
	clog<<"max_IRdrop after SA  :	"<<max_IRdrop<<endl;	
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
	//solve_init();
	// build up two VDD pad sets
	//pad_set_init();
	//for(size_t i=0;i<VDD_set.size();i++)
		//clog<<"i, vdd: "<<i<<" "<<*VDD_set[i]<<endl;	
	//clog<<"after pad set init. "<<endl;
	//solve_LU_core();

	solve_GS();
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
				    !(*it)->ab[1]->is_ground() ){
					continue; // it's a 0v via
				}
				stamp_VDD(A, b, (*it));
			}
			break;
		case PAD:
			break;
		default:
			report_exit("Unknown net type\n");
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
	Node *&rm_pad, Node *&add_pad,  
	size_t &rm_pad_index, size_t iter_move){
	size_t VDD_num = VDD_set.size();
	vector<Node*> nbr_pads;
	size_t nbr_index;

	// the nodes updated in each movement 
	nodesUpdate_move.resize(0);

	//1. randomly pick a Vdd pad to remove
	rm_pad_index = random_gen(0, VDD_num-1); 
	rm_pad = VDD_set[rm_pad_index];

	//Mean_shift_one_move(rm_pad, rm_pad_index, add_pad, add_index);
	//clog<<"rm_pad, add pad: "<<*rm_pad<<" "<<*add_pad<<endl;
	//clog<<"rm_pad: "<<*rm_pad<<endl;
	// It is no more X node
	rm_pad->flag= false;
	rm_pad->value = 0;

	size_t radius = 1;
	add_pad = find_max_IRdrop_candi(rm_pad, radius);
	// find its neiboring pad candidates	
	form_nbr_pads(rm_pad, nbr_pads);	
	if(nbr_pads.size()==0) return;
	//clog<<"nbr_pads.size: "<<nbr_pads.size()<<endl;
	
	// random select a valid nbr pad
	nbr_index = random_gen(0, nbr_pads.size()-1);
	//clog<<"nbr_index: "<<nbr_index<<endl;
	if(add_pad->name == rm_pad->name) 
		add_pad = nbr_pads[nbr_index];
	// choose candi pads by mean shift method
	add_pad->flag = true;
	add_pad->value = VDD;
	solve_GS();
	//clog<<"after solve GS. "<<endl;
	//clog<<"rm_pad, add_pad: "<<*rm_pad<<" "<<*add_pad<<endl; 

	// update pad value and nbr area by iterations
	//update_pad_value(rm_pad, add_pad, nodesUpdate_move, 
			//iter_move);
	
}

void Circuit::one_move_modified(vector<Node*>&nodesUpdate_move,
	double *rhs, Node *&rm_pad, Node *&add_pad,  
	size_t &rm_pad_index, size_t iter_move){
	size_t VDD_num = VDD_set.size();
	vector<Node*> nbr_pads;
	size_t nbr_index;

	// the nodes updated in each movement 
	nodesUpdate_move.resize(0);

	// find its neiboring pad candidates	
	form_nbr_pads(rm_pad, nbr_pads);	
	if(nbr_pads.size()==0) return;
	
	// random select a valid nbr pad
	nbr_index = random_gen(0, nbr_pads.size()-1);
	add_pad = nbr_pads[nbr_index];
	add_pad->flag = true;
	add_pad->value = VDD;

	// update pad value and nbr area by iterations
	update_pad_value(rm_pad, add_pad, nodesUpdate_move, 
			iter_move);
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
					nb->isX()==false){
				nbr_pads.push_back(nb);
				//clog<<"nbr_pads: "<<*nb<<endl;
			}
			else if(rm_pad->name==
				nb->name && na->isX()==false){
				nbr_pads.push_back(na);	
				//clog<<"nbr_pads: "<<*na<<endl;
			}
		}
	}
}

void Circuit::update_pad_value(Node *rm_pad, Node *add_pad, 
	vector<Node*>&nodesUpdate_move, int iter_move){
	size_t iter = 1;
	// eps0 need to be very small
	double eps0=1e-10;
	Node *ref_node;
	ref_node = rm_pad;
		
	int *timestamps;
	timestamps = new int [nodelist.size()-1];
	for(size_t i=0;i<nodelist.size();i++)
		timestamps[i]=0;
	
	size_t MAX_QUEUE_SIZE = nodelist.size()-1;
	
	CircularQueue q(MAX_QUEUE_SIZE);
	double V_improve = 1;
	double V_ref_old = 0;
	double max_diff = 1;
	double diff = 0;
	Node *max_nd;
	// only test 1 iteration
	while(max_diff > eps0){// && iter<LIMIT){
		//V_ref_old = ref_node->value;
		q.reset();
		q.insert(rm_pad);
		q.insert(add_pad);
		rm_pad->flag_visited = iter; 
		add_pad->flag_visited = iter;
		max_diff = 0;
		while(q.isEmpty()==false){
			Node *nd = q.extractFront();
			if(timestamps[nd->rid]!=iter_move){
				nodesUpdate_move.push_back(nd);
				timestamps[nd->rid] = iter_move;
			}
			// update node value with neighbors
			diff = update_node_value(iter, rm_pad, nd);
			if(abs(diff)>max_diff) {
				max_diff =abs(diff);
				max_nd = nd;
			}
			// if there is a change with current nodes'
			// value, queue its neighbors
			// update_queue with nd's neighbors
			if(diff>1e-10)
				update_queue(q, nd, iter);
		}
		//clog<<"iter, max_diff, nd: "<<iter<<" "<<max_diff<<
			//" "<<*max_nd<<endl;
		iter++;
	}
	
	//cout<<"iter for one rm_pad is: "<<iter<<endl;
	delete [] timestamps;
}

void Circuit::update_pad_value_optimize(Node *rm_pad, Node *add_pad, 
	 int iter_move, double *rhs){
	size_t iter = 1;
	double eps0=1e-4;
	double eps1 = 1e-4;
	size_t LIMIT = 500;
	Node *ref_node;
	ref_node = rm_pad;

	rm_pad->value=0;
	add_pad->value=VDD;
	
	int *timestamps;
	timestamps = new int [nodelist.size()-1];
	for(size_t i=0;i<nodelist.size();i++)
		timestamps[i]=0;

	size_t MAX_QUEUE_SIZE = nodelist.size()-1;
	double V_improve = 1;
	double V_ref_old = 0;
	// only test 1 iteration
	while(V_improve > eps0 && iter<LIMIT){
		V_ref_old = ref_node->value;
		//cout<<"vref old: "<<V_ref_old<<endl;
		queue<Node*> q;
		//cout<<"reset q here. "<<endl;
		q.push(rm_pad);
		q.push(add_pad);
		rm_pad->flag_visited = iter; 
		add_pad->flag_visited = iter;
		double diff = 1;
		//size_t count = 0;
		//while(!q.isEmpty() && (V_improve>eps1)){
		while(q.empty()==false){
			Node *nd = q.front();
			if(timestamps[nd->rid]!=iter_move){
				timestamps[nd->rid] = iter_move;
			}
			//count++;
			// update node value with neighbors
			diff = update_node_value(iter, rm_pad, nd);
			if(diff<1e-4){
				//cout<<"iter, nd stop: "<<iter<<" "<<*nd<<endl;
				break;
			}
			// if there is a change with current nodes'
			// value, queue its neighbors
			// update_queue with nd's neighbors
			update_queue_optimize(q, nd, iter);
			//cout<<"qeue is: "<<endl;
			//cout<<q;
			//cout<<endl<<endl;	
		}
		//cout<<"update "<<count<<" nodes. "<<endl;
		V_improve = fabs(V_ref_old - ref_node->value);
		//cout<<"old, new, diff for pad: "<<V_ref_old<<" "
			//<<ref_node->value<<" "<<V_improve<<" "<<
			//*ref_node<<endl;
		iter++;
	}
	//cout<<"iter for one rm_pad is: "<<iter<<endl;
	delete [] timestamps;
}


// rhs is the current vector
double Circuit::update_node_value(int iter, Node *&rm_pad, Node *&nd){
	if(nd->isX()==true) {
		return 0;
	}
	double h = 0.06;
	double omega = 2-h;
	double V_old=0;
	double V_temp = 0;
	double G = 0;
	Net *net;
	Node *nbr, *na, *nb;
	double sum = 0;
	double current = 0;
	net = NULL;
	nbr = NULL; na = NULL; nb = NULL;

	V_old = nd->value;
	
	// update nd->value
	for(int i=0;i<6;i++){
		net = nd->nbr[i];
		if(net ==NULL) continue;
		G = 1.0/net->value;
		na = net->ab[0]; nb = net->ab[1];
		if(nd->name == na->name) nbr = nb;
		else	nbr = na;
		if(!nbr->is_ground()){
			sum += G;
			V_temp += G*nbr->value;
		}
	}
	if(nd->nbr[BOTTOM]== NULL) current = 0;
	else	current = -nd->nbr[BOTTOM]->value;
	V_temp += current;
	V_temp /=sum;
	if(iter ==1 && nd->name==rm_pad->name)
		nd->value  = V_temp;
	else
		nd->value = (1-omega)*nd->value + omega*V_temp;
 	
	double V_improve = fabs(nd->value - V_old);

	return V_improve;
}

void Circuit::update_queue(CircularQueue &q, Node *nd, size_t iter){
	Net * net; Node *nbr;
	Node *na, *nb;
	//cout<<"update queue:head "<<q.queueHead_<<endl;
	//cout<<"center nd: "<<*nd<<endl;
	for(int i=0;i<6;i++){
		net = nd->nbr[i];
		if(net==NULL) continue;
		// find the other node, which is nbr
		na = net->ab[0];
		nb = net->ab[1];
		if(nd->name == na->name)
			nbr = nb;
		else	nbr = na;
		if(!nbr->is_ground()&& nbr->flag_visited != iter){
			q.insert(nbr);
			nbr->flag_visited = iter;
		}
	}
}

void Circuit::update_queue_optimize(queue<Node*> &q, Node *nd, size_t iter){
	Net * net; Node *nbr;
	Node *na, *nb;
	//cout<<"update queue:head "<<q.queueHead_<<endl;
	//cout<<"center nd: "<<*nd<<endl;
	for(int i=0;i<6;i++){
		net = nd->nbr[i];
		if(net==NULL) continue;
		// find the other node, which is nbr
		na = net->ab[0];
		nb = net->ab[1];
		if(nd->name == na->name)
			nbr = nb;
		else	nbr = na;
		if(!nbr->is_ground()&& nbr->flag_visited != iter){
			q.push(nbr);
			nbr->flag_visited = iter;
		}
	}
}

// update cost and total_cost
double Circuit::update_cost(vector<Node*> &nodesUpdate_move, 
		int iter_T, double &change_cost_total,
		double *old_voltages){
	double change_cost = 0;
	Node *nd;
	for (size_t i=0; i<nodesUpdate_move.size();i++){
		nd = nodesUpdate_move[i] ;
		change_cost += penalty(VDD-nd->value, max_IRdrop)-
			penalty(VDD-old_voltages[nd->rid], max_IRdrop);
	}
	if (iter_T==0){
		// accept all the movements at the beginning
		change_cost_total += change_cost;
	}
	return change_cost;
}

// update cost and total_cost
double Circuit::update_cost_modified(vector<Node*> &nodesUpdate_move,
		double *old_voltages){
	double change_cost = 0;
	Node *nd;
	for (size_t i=0; i<nodesUpdate_move.size();i++){
		nd = nodesUpdate_move[i] ;
		change_cost += penalty(VDD-nd->value, max_IRdrop)-
			penalty(VDD-old_voltages[nd->rid], max_IRdrop);
	}
	return change_cost;
}

void Circuit::accept_move(vector<Node*>&nodesUpdate_move, 
		double *old_voltages, size_t rm_index, 
		Node *add_pad){
	Node *nd;
	//size_t size = nodesUpdate_move.size()-1;
	for (size_t i=0; i<nodelist.size()-1;i++){
		nd = nodelist[i];//nodesUpdate_move[i];
		//cout<<"i, rid, nd: "<<i<<" "<<nd->rid<<" "<<*nd<<endl;
		//if (old_voltages[nd->rid] != nd->value) {
			//if (nd->critical==true){
				// remove the node from critical 
				// set O(logN)
				//CriticalNodes.erase(nd) ; 
				//nd->critical = false ; 
			//}
			//old_voltages[nd->rid] = nd->value;
			old_voltages[i] = nd->value; 
			//if ((VDD - nd->value) > th_IRdrop) {
				// insert it into the critical set 
				// O(logN)
				//CriticalNodes.insert(nd) ; 
				//nd->critical = true ;
			//}
		//}
	}
	//keep the change of Vdd
	VDD_set[rm_index]=add_pad;  
}

void Circuit::reject_move(vector<Node*>&nodesUpdate_move, 
	Node *rm_pad, Node *add_pad, double *old_voltages){
	Node *nd;
	for (size_t i=0; i<nodelist.size()-1;i++){//nodesUpdate_move.size();i++){
		nd = nodelist[i];//nodesUpdate_move[i];
		nd->value = old_voltages[i];//nd->rid];
	}
	rm_pad->flag = true; // rm_pad isX again
	rm_pad->value = VDD;
	add_pad->flag = false; //add_pad is not X again
}

void Circuit::recompute_worst_IRdrop(){
	Node *nd;
	set<Node*, VoltageLessThan>::iterator it;
	if (!CriticalNodes.empty()) {
		max_IRdrop=0;
		for(it = CriticalNodes.begin();it!=
			CriticalNodes.end();it++){
			nd = *it ;
			if((VDD-nd->value)>max_IRdrop)
				max_IRdrop = VDD-nd->value;
		}
	}else {
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
    return min + int( ((max-min) +1) * (rand() /(RAND_MAX + 1.0))); 
}

double Circuit::optimize_pad_assign(){
	Node *rm_pad, *add_pad;
	double cost = 0;
	size_t Movement = 10;
	double final_cost = 0;
	rm_pad = NULL;
	add_pad = NULL;
	size_t min_index = 0;

	vector<Node *> nodesUpdate_move;
		
	for(size_t iter=1; iter<Movement; iter++){	
		nodesUpdate_move.resize(0);
		// find pad located in mimum IR drop area
		rm_pad = find_min_IRdrop_pad(min_index);
		rm_pad->flag = false; // rm_pad is not VDD pad now
		// find candidate pad located in maximum IR drop area
		add_pad = find_max_IRdrop_candi();
		add_pad->flag = true; // add_pad is VDD pad now
		add_pad->value = VDD;
		rm_pad->value = 0;
		//cout<<"iter rm_pad add_pad: "<<iter<<" "<<
			//*rm_pad<<" "<<*add_pad<<endl;
		//cout<<"rm_pad "<<*rm_pad<<endl;
		//cout<<"add_pad "<<*add_pad<<endl;

		update_pad_value(rm_pad, add_pad, nodesUpdate_move, 
				iter);

		VDD_set[min_index] = add_pad;
	}
	locate_maxIRdrop();
	locate_thIRdrop();
	final_cost = SACost();
	//clog<<"final_cost: "<<final_cost<<endl;	
	nodesUpdate_move.clear();
	return final_cost;
}

double Circuit::optimize_pad_assign_new(){
	Node *rm_pad, *add_pad;
	double cost = 0;
	size_t Movement = 20;
	double final_cost = 0;
	rm_pad = NULL;
	add_pad = NULL;
	size_t min_index = 0;
	size_t iter = 0;
	Node *prev_rm_pad, *prev_add_pad;
	rm_pad = NULL; add_pad = NULL;
	prev_rm_pad = NULL;
	prev_add_pad = NULL;

	vector<Node *> nodesUpdate_move;
	while(1){
		clear_candi_visit();
	//for(size_t iter=1; iter<Movement; iter++){	
		nodesUpdate_move.resize(0);
		// find pad located in mimum IR drop area
		rm_pad = find_min_IRdrop_pad(min_index);
		add_pad = find_max_IRdrop_candi();
		rm_pad->flag_visited ++;
		
		// local bumping shows up, break
		if(prev_add_pad!=NULL && prev_rm_pad !=NULL && 
			(rm_pad->name == prev_add_pad->name && 
		   	add_pad->name == prev_rm_pad->name))
			// for untouched pads flag_visited != 0
			break;
		//clog<<"rm_pad, add_pad: "<<*rm_pad<<" "<<*add_pad<<endl;	
		prev_rm_pad = rm_pad;
		prev_add_pad = add_pad;
		
		add_pad->flag = true; // add_pad is VDD pad now
		add_pad->value = VDD;
		rm_pad->flag = false;
		rm_pad->value = 0;
		//cout<<"iter rm_pad add_pad: "<<iter<<" "<<
			//*rm_pad<<" "<<*add_pad<<endl;
		//clog<<"rm_pad "<<*rm_pad<<endl;
		//clog<<"add_pad "<<*add_pad<<endl;

		update_pad_value(rm_pad, add_pad, nodesUpdate_move, 
				iter);

		VDD_set[min_index] = add_pad;
	}
	locate_maxIRdrop();
	locate_thIRdrop();
	final_cost = SACost();
	//clog<<"final_cost after initial reallocation: "<<final_cost<<endl;	
	nodesUpdate_move.clear();

	return final_cost;
}


// loop the VDD pads, to find out one that has minimum 
// IR drop area
Node* Circuit::find_min_IRdrop_pad(size_t & min_index){
	Node *nd;
	double IRdrop;
	double minIRdrop = 0;
	Node *min_pad;
	nd = NULL; min_pad = NULL;
	// judge by sum of nbr candi and pad IRdrops
	for(size_t i=0;i<VDD_set.size();i++){
		nd = VDD_set[i];
		IRdrop = area_IRdrop(nd);
		//cout<<"nd, IRdrop: "<<*nd<<" "<<IRdrop<<endl;
		if(i==0){
			minIRdrop = IRdrop;
			min_pad = nd;
			min_index = 0;
		}else if(IRdrop < minIRdrop){
			minIRdrop = IRdrop;
			min_pad = nd;
			min_index = i;	
		}
	}
	return min_pad;
}

// find area IRdrop around node nd
double Circuit::area_IRdrop(Node *nd){
	double sum = 0;
	Node *nbr, *na, *nb;
	Net *net;
	size_t count = 0;
	
	for(int i=0;i<4;i++){
		net = nd->nbr_pad[i];
		if(net== NULL) continue;
		na = net->ab[0];
		nb = net->ab[1];
		if(na->name == nd->name)
			nbr = nb;
		else nbr = na;
		count++;

		sum += VDD - nbr->value;
	}
	return sum/count;
}

Node* Circuit::find_max_IRdrop_candi(){
	Node *nd;
	double IRdrop;
	double maxIRdrop = 0;
	Node *max_pad;
	nd = NULL; max_pad = NULL;
	// judge by sum of nbr candi and pad IRdrops
	for(size_t i=0;i<VDD_candi_set.size();i++){
		nd = VDD_candi_set[i];
		IRdrop = VDD - nd->value;
		if(IRdrop > maxIRdrop){
			maxIRdrop = IRdrop;
			max_pad = nd;
		}	
	}
	return max_pad;
}

void Circuit::rebuild_voltage_nets(){
	int type = VOLTAGE;
	net_set[type].clear();
	for(size_t i=0;i<VDD_set.size();i++){
		Net *net = new Net(VOLTAGE, VDD, VDD_set[i], 
			nodelist[nodelist.size()-1]);		
		net_set[type].push_back(net);	
	}
}

void Circuit::solve_GS(double *rhs){
	double max_diff = 1;
	int iter = 0;
	double omega= 2-0.06;
	Node *max_nd;
	//clog<<"nodelist.size: "<<nodelist.size()-1<<endl;
	while(max_diff >1e-8){// && iter < 500){
		max_diff = 0;
		for(size_t i=0;i<nodelist.size()-1;i++){
			Node *nd = nodelist[i];
			
			if(nd->isX()==true) continue;
			double V_old=0;
			double V_temp = 0;
			double G = 0;
			Net *net;
			Node *nbr, *na, *nb;
			double sum = 0;
			net = NULL;
			nbr = NULL; na = NULL; nb = NULL;

			V_old = nd->value;

			// update nd->value
			for(int i=0;i<6;i++){
				net = nd->nbr[i];
				if(net ==NULL) continue;
				G = 1.0/net->value;
				na = net->ab[0]; nb = net->ab[1];
				if(nd->name == na->name) nbr = nb;
				else	nbr = na;
				if(!nbr->is_ground()){
					sum += G;
					V_temp += G*nbr->value;
				}
			}

			V_temp += rhs[nd->rid];
			V_temp /=sum;
				nd->value = (1-omega)*nd->value + omega * V_temp;

			double diff = fabs(nd->value - V_old);
			if(diff > max_diff) {
				max_diff = diff;
				max_nd = nd;
			}

			//cout<<"nd, v_old, v_temp, diff: "<<nd->name<<" "<<V_old<<" "<<V_temp<<" "<<diff<<endl;
		}
		//cout<<"iter, max_diff: "<<iter<<" "<<max_diff<<endl;
		//cout<<"update "<<count<<" nodes. "<<endl;
		//V_improve = fabs(V_ref_old - ref_node->value);
		iter++;
	}

}

void Circuit::clear_candi_visit(){
	Node *nd;
	for(size_t i=0;i<VDD_candi_set.size();i++){
		nd = VDD_candi_set[i];
		nd->flag_visited = 0;
	}
}

// produce one set of VDD Pad 
void Circuit::random_init_iter(){
	//for(size_t i=0;i<VDD_candi_set.size();i++)
		//if(VDD_candi_set[i]->isX()==true)
			//clog<<"VDD_candi_set. flag: "<<*VDD_candi_set[i]<<endl;
	// clear all X node
	for(size_t i=0;i<VDD_set.size();i++){
		VDD_set[i]->flag = false;
	}
	size_t VDD_num = VDD_set.size();
	//VDD_set.clear();
	size_t index=0;
	for(size_t i=0;i<VDD_num;i++){
		while(1){
			index = random_gen(0,VDD_candi_set.size()-1);
			//clog<<"index, candi: "<<index<<" "<<VDD_candi_set[index]->isX()<<endl;
			if(VDD_candi_set[index]->isX()!= true){
				//clog<<"push this node. "<<*VDD_candi_set[index]<<endl;
				//VDD_set.push_back(VDD_candi_set[index]);
				VDD_set[i]=VDD_candi_set[index];
				VDD_set[i]->flag = true;
				VDD_set[i]->value = VDD;
				//VDD_candi_set[index]->flag= true;
				break;
			}	
		}
	}

}

void Circuit::RANSAC_init(){
	vector<Node*> best_VDD_set;
	best_VDD_set = VDD_set;	
	double best_maxIRdrop = max_IRdrop;
	// randomly select p patterns for init distribution
	for(size_t i=0;i<50;i++){
		// produce one pad assignment
		random_init_iter();

		solve_GS();
		// then solve Lu after rebuild matrix
		//rebuild_voltage_nets();
		//solve_LU_core();
		locate_maxIRdrop();
		double current_maxIRdrop = max_IRdrop;
		if(abs(current_maxIRdrop)<abs(best_maxIRdrop)){
			best_VDD_set = VDD_set;
			// set the VDD_set and VDD_candi_set to be current
			//clog<<"i, accept, current, best: "<<i<<" "<<current_maxIRdrop<<" "<<best_maxIRdrop<<endl;

			best_maxIRdrop = current_maxIRdrop;
		}
	}

	// clear all X node
	for(size_t i=0;i<VDD_set.size();i++){
		VDD_set[i]->flag = false;
	}
	VDD_set = best_VDD_set;
	for(size_t i=0;i<VDD_set.size();i++){
		//clog<<"VDD best, i, node' "<<i<<" "<<*VDD_set[i]<<endl;
		VDD_set[i]->flag = true;
		VDD_set[i]->value = VDD;
	}
	solve_GS();
	//rebuild_voltage_nets();
	//solve_LU_core();
	best_VDD_set.clear();
}

void Circuit::opti_SA(){
	double cost = 0;
	// iteration may broke the single descending characteristics	
	for(size_t i=0;i<1;i++){
		cost = optimize_pad_assign_new();
		
		locate_maxIRdrop();
		clog<<endl<<" max_IRdrop: "<<max_IRdrop<<endl;
		SA(0.01);
		locate_maxIRdrop();
		//locate_maxIRdrop();
		clog<<"max_IRdrop: "<<max_IRdrop<<endl;	
	}
}

// perform mean shift like movement of pads
void Circuit::Mean_shift_move(){
	// define number of nodes in the window for centroid
	int N = 1;
	Node *add_pad;
	size_t add_candi_index;
	// N iterations
	for(size_t i=1;i<=N;i++){
		for(size_t j=0;j<VDD_set.size();j++){
			Node *rm_pad = VDD_set[j];	
			Mean_shift_one_move(rm_pad, j, add_pad, 
				add_candi_index);
			
			rebuild_voltage_nets();
			solve_LU_core();
			// then update node values
		}
		// update node values after a move
	}
}

// for a pad, find its new location in candi set
void Circuit::Mean_shift_one_move(Node *center, size_t center_index, 
	Node *&new_center, size_t &new_candi_index){
	// define number of nodes in the window for centroid
	int num_LIMIT = 500;
	Node *nd, *na, *nb, *nbr;
	double weight = 0;
	double weight_c;
	long mean_x, mean_y, mean_z;
	double diff_x, diff_y, diff_z;
	size_t min_index;
	Net *net;
	CircularQueue q(nodelist.size()-1);
	
	for(size_t j=0;j<nodelist.size()-1;j++)
		nodelist[j]->flag_visited = 0;
	int num_nodes = 0;
	q.reset();
	q.insert(center);
	center->flag_visited = 1;
	//cout<<endl<<"center node: "<<*center<<endl;
	weight = VDD/center->value;
	double sum_x = 0;
	double sum_y = 0;
	double sum_z = 0;
	num_nodes ++;

	while(q.isEmpty()==false && num_nodes <= num_LIMIT){
		nd = q.extractFront();
		//clog<<"nd: "<<*nd<<endl;
		for(size_t j=0;j<6;j++){
			net = nd->nbr[j];
			if(net == NULL) continue;
			na = net->ab[0]; nb = net->ab[1];
			if(nd->name == na->name) nbr = nb;
			else	nbr = na;
			if(!nbr->is_ground()&& nbr->flag_visited ==0){
				//cout<<"nbr: "<<*nbr<<endl;
				weight_c = VDD/nbr->value;
				diff_x = (nbr->pt.x - center->pt.x);
				diff_y = (nbr->pt.y - center->pt.y);
				diff_z = (nbr->pt.z - center->pt.z);
				//cout<<"diff_x, y, z: "<<diff_x<<" "<<diff_y<<" "<<diff_z<<endl;
				sum_x += diff_x*weight_c;
				sum_y += diff_y*weight_c;
				sum_z += diff_z*weight_c;

				//cout<<"sum_x, y, z: "<<sum_x<<" "<<sum_y<<" "<<sum_z<<endl;
				q.insert(nbr);
				nbr->flag_visited ++;
				num_nodes ++;
			}
		}
	}
	//clog<<"num_nodes: "<<num_nodes<<endl;	
	// after getting the 100 nodes sum, compute mean
	mean_x = (long)(center->pt.x + sum_x / num_nodes);
	mean_y = (long)(center->pt.y + sum_y / num_nodes);
	mean_z = (long)(center->pt.z + sum_z / num_nodes);
	//cout<<center->pt.x+sum_x/num_nodes<<" "<<center->pt.y + 
		//sum_y / num_nodes<<" "<<center->pt.z + sum_z / num_nodes<<endl;
	//cout<<"new center: "<<mean_x<<" "<<mean_y<<" "<<mean_z<<endl;
	// map this new center into one of the candidate location
	double dist, min_dist;
	size_t count = 0;
	// scan all candi location to find the closest one
	for(size_t k=0;k<VDD_candi_set.size();k++){
		Node * candi = VDD_candi_set[k];
		if(candi->isX()== true) continue;
		count ++;
		//cout<<endl<<"k, candi: "<<k<<" "<<*candi<<endl;
		//cout<<"candi->pt.x, mean_x: "<<candi->pt.x<<" "<<
		//mean_x<<" "<<candi->pt.x-mean_x<<endl;
		diff_x = fabs(candi->pt.x-mean_x);
		diff_y = fabs(candi->pt.y-mean_y);
		diff_z = fabs(candi->pt.z-mean_z);
		dist = sqrt(diff_x*diff_x + diff_y*diff_y + 
				diff_z*diff_z); 
		//cout<<"min_dist, dist: "<<min_dist<<" "<<dist<<endl;
		if(count==1){
			min_dist = dist;
			min_index = k;
		}
		else if(dist < min_dist){
			min_dist = dist;
			min_index = k;
		}
	}

	//cout<<"min_dist, index: "<<min_dist<<" "<<min_index<<endl;

	//cout<<"corresponding candi: "<<*VDD_candi_set[min_index]<<endl;
	// update node flag in VDD_set
	//VDD_set[center_index]->flag = false;
	//VDD_set[center_index]->value = 0;
	//VDD_set[center_index] = VDD_candi_set[min_index];
	//VDD_set[center_index]->flag = true;
	//VDD_set[center_index]->value = VDD;
	new_center = VDD_candi_set[min_index];
	new_candi_index = min_index;
	q.reset();
}

// find the maximum candi locations within radius of rm_pad
Node* Circuit::find_max_IRdrop_candi(Node *rm_pad, size_t Limit){
	Node *nd;
	double IRdrop;
	double maxIRdrop = 0;
	Node *max_pad;
	Node *na, *nb, *nbr;
	nd = NULL; max_pad = NULL;
	// store the original flag_visited for nodes
	//vector<size_t> VDD_candi_flag_visited_temp;
	//VDD_candi_flag_visited_temp.resize(VDD_candi_set.size());
	//for(size_t i=0;i<VDD_candi_set.size();i++){
		//VDD_candi_flag_visited_temp[i] = 
			//VDD_candi_set[i]->flag_visited;
		//VDD_candi_set[i]->flag_visited = 0;
	//}
	//size_t count = 0;
	nd = rm_pad;
	//rm_pad->flag_visited ++;
	//clog<<endl<<"rm_pad: "<<*rm_pad<<endl;
	//while(nd != NULL && count < Limit){
		//clog<<"center node: "<<*nd<<endl;
		for(size_t i=0;i<4;i++){
			Net *net = nd->nbr_pad[i];
			if(net== NULL) continue;
			na = net->ab[0];
			nb = net->ab[1];
			if(na->name == nd->name && !nb->is_ground())
				nbr = nb;
			else if(nb->name == nd->name && !na->is_ground())
				nbr = na;
			//if(nbr->flag_visited != 0) continue;
			//clog<<"effective nbr: "<<*nbr<<endl;
			//count ++;
			if(VDD-nbr->value > maxIRdrop){
				maxIRdrop = VDD-nbr->value;
				max_pad = nbr;
				//nd = nbr;
				//nbr->flag_visited ++;				
			}
		}
	//}	
	// assign flag_visited back
	//for(size_t i=0;i<VDD_candi_set.size();i++)
		//VDD_candi_set[i]->flag_visited = 
		  //VDD_candi_flag_visited_temp[i];
	//VDD_candi_flag_visited_temp.clear();
	return max_pad;
}

void Circuit::solve_GS(){
	double max_diff = 1;
	int iter = 0;
	double omega= 2-0.06;
	Node *max_nd;
	//clog<<"nodelist.size: "<<nodelist.size()-1<<endl;
	while(max_diff >1e-8){// && iter < 500){
		max_diff = 0;
		for(size_t i=0;i<nodelist.size()-1;i++){
			Node *nd = nodelist[i];
			
			if(nd->isX()==true) continue;
			double V_old=0;
			double V_temp = 0;
			double G = 0;
			Net *net;
			Node *nbr, *na, *nb;
			double sum = 0;
			net = NULL;
			nbr = NULL; na = NULL; nb = NULL;

			V_old = nd->value;

			// update nd->value
			for(int i=0;i<6;i++){
				net = nd->nbr[i];
				if(net ==NULL) continue;
				G = 1.0/net->value;
				na = net->ab[0]; nb = net->ab[1];
				if(nd->name == na->name) nbr = nb;
				else	nbr = na;
				if(!nbr->is_ground()){
					sum += G;
					V_temp += G*nbr->value;
				}
			}

			V_temp += -nd->nbr[BOTTOM]->value;//rhs[nd->rid];
			V_temp /=sum;
				nd->value = (1-omega)*nd->value + omega * V_temp;

			double diff = fabs(nd->value - V_old);
			if(diff > max_diff) {
				max_diff = diff;
				max_nd = nd;
			}

			//cout<<"nd, v_old, v_temp, diff: "<<nd->name<<" "<<V_old<<" "<<V_temp<<" "<<diff<<endl;
		}
		//clog<<"iter, max_diff: "<<iter<<" "<<max_diff<<endl;
		//cout<<"update "<<count<<" nodes. "<<endl;
		//V_improve = fabs(V_ref_old - ref_node->value);
		iter++;
	}

}
