#include <fstream>
#include <string>
#include <vector>
#include "mg_circuit.h"
using namespace std;

MG_Circuit::MG_Circuit(string name){
	name=name;
	LEVEL=0;		
	mg_ckt.clear();
}

MG_Circuit::~MG_Circuit(){
	mg_ckt.clear();
}

void build_mg_ckt(Circuit *ckt){
}

void solve_mg_ckt(Circuit *ckt){
} 
