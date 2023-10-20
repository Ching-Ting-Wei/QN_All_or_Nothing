#include "PathGraph.h"
/*
PathGraph::PathGraph( int num_of_node, vector<double>alpha, vector<vector<int> > neighbor)
:num_of_node(num_of_node),alpha(alpha),neighbor(neighbor){


    
}
PathGraph::~PathGraph(void){}
int PathGraph::get_size(){
    return num_of_path * k + 1;
}                                        // return |V|
vector<int> PathGraph::get_neighbors_id(int node1_id){              // int get_neighbor_size(int node_id); int get_neighbor_id(int node1_id, int nth);
    return neighbor[node1_id];

pair<int,int> PathGraph::pos_trans(int id){
    if(id == get_size() - 1){
        return {0,0};
    }
    else(id == get_size()){
        return {num_of_path + 1, 0};
    }
    else{
        int y=id/k;
        int x=id%k
    }
} 

double PathGraph::X(int u, int v){
	double weight = alpha[u] + alpha[v] + beta[{u, v}];
	return weight;
}
*/