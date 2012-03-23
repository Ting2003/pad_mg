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
		ss << i;
		name = ss.str();
		mg_ckt[i] = new Circuit(name);
		if(i==0){
			clog<<"start build mg_ckt[0]. "<<endl;
			mg_ckt[i] = build_one_layer_circuit(ckt);
		}
		else
			mg_ckt[i] = build_one_layer_circuit(mg_ckt[i]);
	}
}

// build coarser grid from fine one
Circuit * MG_Circuit::build_one_layer_circuit(Circuit *ckt){
	Circuit* coarse_ckt;
	coarse_ckt = new Circuit();
	coarse_ckt->nodelist.clear();
	// build nodelist of coarse_ckt
	Point pt_c, prev_pt; 
	int count_x=0;
	int count_y=0;
	// extract ground node exclusively
	for(size_t i=0;i<ckt->nodelist.size()-1;i++){
		Node *nd = ckt->nodelist[i];
		//clog<<endl<<"ckt_nd: "<<*nd<<endl;
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
		//clog<<"count_x and y: "<<count_x<<" "<<count_y<<endl;
		// keep this node in coarse grid
		if(count_y % 2==0 && count_x % 2 ==0){			
			pt_c.x = nd->pt.x / 2;
			pt_c.y = nd->pt.y / 2;
			// suppose there is only 1 layer in z direction
			pt_c.z = nd->pt.z ;

			// add this node into coarse nodelist
			Node *nd_c = new Node(nd->name, pt_c, false, 0.0);
			//clog<<"nd_c: "<<*nd_c<<endl;
			coarse_ckt->nodelist.push_back(nd_c);
			//clog<<"push back nd_c. "<<endl;
			count_y = 0;
			count_x = 0;	
		}
	}
	//clog<<"before ground node. "<<endl;
	// handle ground node
	Node *nd = ckt->nodelist[ckt->nodelist.size()-1];
	Node *nd_c = new Node(*nd);
	coarse_ckt->nodelist.push_back(nd_c);
	//for(size_t i=0;i< coarse_ckt->nodelist.size();i++){
		//cout<<"i, coarse_ckt_nodes: "<<i<<" "<<
		  //*coarse_ckt->nodelist[i]<<" "<<coarse_ckt->nodelist[i]->pt<<endl;	
	//}
	return coarse_ckt;
}

//void MG_Circuit::solve_mg_ckt(Circuit *ckt){
//} 
