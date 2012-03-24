#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include "mg_circuit.h"
using namespace std;

MG_Circuit::MG_Circuit(){
	LEVEL=0;		
	mg_ckt.clear();
}

MG_Circuit::~MG_Circuit(){
	mg_ckt.clear();
}

// build Level layers of circuits
void MG_Circuit::build_mg_ckt(Circuit *ckt, int layer){
	LEVEL = layer;
	// allocate layer of circuits
	mg_ckt.resize(LEVEL);
	stringstream ss;
	string name;
	for(int i=0;i<LEVEL;i++){
		//ss << i;
		//name = ss.str();
		//mg_ckt[i] = new Circuit(name);
		if(i==0){
			//clog<<"start build mg_ckt[0]. "<<endl;
			build_one_layer_circuit(ckt, i);
		}
		else
			build_one_layer_circuit(mg_ckt[i-1], i);
	}
}

// build coarser grid from fine one
Circuit * MG_Circuit::build_one_layer_circuit_nodelist(Circuit *ckt){
	Circuit *coarse_ckt;
	coarse_ckt = new Circuit();
	coarse_ckt->nodelist.clear();
	// build nodelist of coarse_ckt
	Point pt_c, prev_pt; 
	int count_x=0;
	int count_y=0;
	// extract ground node exclusively
	for(size_t i=0;i<ckt->nodelist.size()-1;i++){
		Node *nd = ckt->nodelist[i];
		//clog<<endl<<"center node: "<<*nd<<endl;
		if(i==0){
			prev_pt = nd->pt;
		}
		else{
			if(nd->pt.x == prev_pt.x && nd->pt.z == prev_pt.z)
				count_y ++;
			if(nd->pt.y == prev_pt.y && nd->pt.z == prev_pt.z)
				count_x++;
		}
		if(count_y ==2)
			prev_pt = nd->pt;
		// keep this node in coarse grid
		if(count_y % 2==0 && count_x % 2 ==0){			
			pt_c.x = nd->pt.x / 2;
			pt_c.y = nd->pt.y / 2;
			// suppose there is only 1 layer in z direction
			pt_c.z = nd->pt.z ;

			// add this node into coarse nodelist
			Node *nd_c = new Node(nd->name, pt_c, false, 0.0);
			coarse_ckt->add_node(nd_c);
			count_y = 0;
			count_x = 0;	
		}
	}
	// handle ground node
	Node *nd = ckt->nodelist[ckt->nodelist.size()-1];
	Node *nd_c = new Node(*nd);	
	// add G into nodelist, and build up map
	coarse_ckt->add_node(nd_c);

	
	//for(size_t i=0;i< coarse_ckt->nodelist.size();i++){
		//cout<<"i, coarse_ckt_nodes: "<<i<<" "<<
		  //*coarse_ckt->nodelist[i]<<" "<<coarse_ckt->nodelist[i]->pt<<endl;	
	//}
	return coarse_ckt;
}

// ckt is fine grid
void MG_Circuit::build_one_layer_circuit(Circuit *ckt, int level){
	mg_ckt[level] = build_one_layer_circuit_nodelist(ckt);
	Node *nd, *nd_c;
	//clog<<mg_ckt[level]->nodelist.size()<<endl;
	//for(size_t i=0;i< mg_ckt[level]->nodelist.size()-1;i++){
		//cout<<"i, coarse_ckt_nodes: "<<i<<" "<<
		  //*mg_ckt[level]->nodelist[i]<<endl;	
	//}
	for(size_t i=0;i<mg_ckt[level]->nodelist.size()-1;i++){
		nd_c = mg_ckt[level]->nodelist[i];
		// find node in finer grid
		nd = ckt->get_node(nd_c->name);
		// build the nbr nets
		set_nbr_nets(nd, nd_c, ckt, mg_ckt[level]);
		set_VDD_pads(ckt, mg_ckt[LEVEL]);
	}
}
// for each node in coarse grid, build its neighboring nets 
// also build current nets on coarse grid
void MG_Circuit:: set_nbr_nets(Node *nd, Node *&nd_c, Circuit *ckt,
	Circuit *&coarse_ckt){
	Net *nbr_net;
	Node *center, *nd_a, *nd_b;
	Net *nbr_net_x, *nbr_net_y, *nbr_net_c;
	Node *nbr_x, *nbr_y;
	Node *nbr_xc, *nbr_yc;
	Net *coarse_net_x;
	Net *coarse_net_y;
	Net *coarse_net_current;
	double value_x, value_y, current;
	value_x = 0;
	value_y = 0;
	current = 0;
	nbr_net = NULL; center = NULL; nd_a = NULL; nd_b = NULL;
	nbr_net_x = NULL; nbr_net_y = NULL; nbr_net_c = NULL;
	coarse_net_x = NULL; coarse_net_y = NULL;
	coarse_net_current = NULL;
	// store the 4 interest node;
	Node *nd_temp[4];
	for(int i=0;i<4;i++)
		nd_temp[i] = NULL;
	nd_temp[0] = nd;
	for(int i=1;i<4;i++){
		if(i<3) center = nd;
		else	center = nd_temp[i-1];
		if(center == NULL) continue;
		//clog<<"i, center: "<<i<<" "<<*center<<endl;
		if(i==1 || i==3)
			nbr_net = center->nbr[EAST];
		else	nbr_net = center->nbr[NORTH];
		if(nbr_net == NULL) continue;
		//clog<<"nbr net: "<<*nbr_net<<endl;
		nd_a = nbr_net->ab[0];
		nd_b = nbr_net->ab[1];
		//clog<<"center, nd_a, nd_b: "<<*center<<" "<<*nd_a<<" "<<*nd_b<<endl;
		if(nd_a->name == center->name) 
			nd_temp[i] = nd_b;
		else
			nd_temp[i] = nd_a;
		//clog<<"i, nd_temp[i]: "<<i<<" "<<*nd_temp[i]<<endl;
	}
	
	// then form the nbr resistor and current nets
	for(int i=0;i<4;i++){
		if(nd_temp[i]==NULL) continue;
		center = nd_temp[i];
		nbr_net_x = center->nbr[EAST];
		nbr_net_y = center->nbr[NORTH];
		nbr_net_c = center->nbr[BOTTOM];
		if(nbr_net_x != NULL)
			value_x += nbr_net_x->value;
		if(nbr_net_y != NULL)
			value_y += nbr_net_y->value;
		if(nbr_net_c != NULL)
			current += nbr_net_c->value;	
	}
	

	nbr_x = NULL; nbr_y = NULL;
	nbr_xc = NULL; nbr_yc = NULL;
	if(value_x != 0 && nd_temp[1]!=NULL){
		nbr_net = nd_temp[1]->nbr[EAST];
		if(nbr_net !=NULL){
			nd_a = nbr_net->ab[0];
			nd_b = nbr_net->ab[1];
			if(nd_a->name == nd_temp[1]->name)
				nbr_x = nd_b;
			else 	nbr_x = nd_a;
		}
	}
	if(value_y != 0 && nd_temp[2]!=NULL){
		nbr_net = nd_temp[2]->nbr[NORTH];
		if(nbr_net !=NULL){
			nd_a = nbr_net->ab[0];
			nd_b = nbr_net->ab[1];
			if(nd_a->name == nd_temp[2]->name)
				nbr_y = nd_b;
			else 	nbr_y = nd_a;
		}
	}
	
	// nbr_x and nbr_y are nodes in fine grid
	// nbr_xc and nbr_yc are nodes in coarse grid
	if(nbr_x != NULL){
		nbr_xc = coarse_ckt->get_node(nbr_x->name);
		if(nbr_xc != NULL){
			coarse_net_x = new Net(RESISTOR, value_x, 
				nd_c, nbr_xc);
		}
	}
	if(nbr_y != NULL){
		nbr_yc = coarse_ckt->get_node(nbr_y->name);
		if(nbr_yc != NULL){
			coarse_net_y = new Net(RESISTOR, value_y, 
				nd_c, nbr_yc);
		}
	}
	if(current != 0){
		Node *Ground = coarse_ckt->nodelist[
			coarse_ckt->nodelist.size()-1];
		coarse_net_current = new Net(CURRENT, current, nd_c, Ground);
	}
	// assign these nets into nodes->nbr
	// nbr_x, nbr_y  is nd's nbrs in coarse grid
	nd_c->nbr[EAST] = coarse_net_x;
	nd_c->nbr[NORTH] = coarse_net_y;
	nd_c->nbr[BOTTOM] = coarse_net_current;
	
	if(nbr_xc != NULL)	
		nbr_xc->nbr[WEST] = nd_c->nbr[EAST];//coarse_net_x;
	if(nbr_yc != NULL)
		nbr_yc->nbr[SOUTH] = nd_c->nbr[NORTH];//coarse_net_y;
	
	//if(nd_c->nbr[EAST] != NULL)
		//clog<<"new net_x: "<<*nd_c->nbr[EAST]<<endl;
	//if(nd_c->nbr[NORTH] != NULL)
		//clog<<"new net_y: "<<*nd_c->nbr[NORTH]<<endl;
	//if(nd_c->nbr[BOTTOM]!=NULL)
		//clog<<"new current: "<<*nd_c->nbr[BOTTOM]<<endl;
}

// set VDD_set and VDD_candi_set
void MG_Circuit::set_VDD_pads(Circuit *ckt, Circuit *&coarse_ckt){
	coarse_ckt->VDD_set.resize(ckt->VDD_set.size());
	coarse_ckt->VDD_candi_set.resize(ckt->VDD_candi_set.size());
	Node *nd, *nd_c;
	for(size_t i=0;i<ckt->VDD_set.size();i++){
		nd = ckt->VDD_set[i];
		nd_c = coarse_ckt->get_node(nd->name);
		coarse_ckt->VDD_set[i] = nd_c;
	}
	for(size_t i=0;i<ckt->VDD_candi_set.size();i++){
		nd = ckt->VDD_candi_set[i];
		nd_c = coarse_ckt->get_node(nd->name);
		coarse_ckt->VDD_candi_set[i] = nd_c;
	}
}
//void MG_Circuit::solve_mg_ckt(Circuit *ckt){
//} 
