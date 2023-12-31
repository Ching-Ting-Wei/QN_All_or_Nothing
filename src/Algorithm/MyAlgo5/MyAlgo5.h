#ifndef __MyAlgo5_H
#define __MyAlgo5_H

#include <iostream>
#include <algorithm>
#include <omp.h>
#include <queue>
#include <limits>
#include <string>
#include <cmath>
#include <iomanip>
#include "../AlgorithmBase/AlgorithmBase.h"
#include "../../Request/WholeRequest.h"
#include "../../Network/Graph/Graph.h"
#include "./PathGraph.h"
#include "../../config.h"

using namespace std;

class MyAlgo5:public AlgorithmBase {

private:
    vector<map<pair<int,int>, double>> Y;
    vector<double> alpha;
    map<vector<int>, double> x_i_p;
    vector<map<vector<int>, double>> x_i_s;
    map<pair<int,int>, double> beta;
    vector<vector<vector<int>>> all_source_target_path;
    vector<vector<vector<int>>> all_given_path;
    int change_edge_num = 0; 
    int diff_num = 0;
    int qubit_num = 3;
    double epsilon = 0.2;
    int path_num ;
    double obj = 0;
    double delta;
    double M;  
    vector<double> tau;               
    double X(int u, int v, int i);
public: 
    map<int, int> num_of_path_count;
    map<int, int> path_length_encode;
    map<int, int> path_length_cnt;
    vector<vector<int>> allPathsSourceTarget(int src, int dst);
    vector<int> SepDijkstra(int dst, int req_no, vector<vector<double>> &path_graph_X);
    vector<int> Dijkstra(vector<vector<int>>&copy_graph, int src, int dst, vector<int>&get_delete);
    vector<int> separation_oracle(int req_no, double &req_Us, vector<vector<vector<double>>> &path_graph_X, vector<vector<vector<double>>> &path_graph_Y);
    vector<map<vector<int>, int>> Greedy_rounding();
    void path_assignment();
    void calculate();
    void entangle();
    void swap();
    void send();
    void dfs(int s, int t, vector<vector<int>> &ans, vector<int> &path, vector<bool> &visited);    
    void next_time_slot();
    void yen(int src,int dst,int K,int req_no);
    void find_bottleneck(vector<int>, int req_no);
    void initialize();
    void check_enough(vector<map<vector<int>, int>> &path);
    void readd(vector<map<vector<int>, int>> &path,vector<int> &over_memory,map<vector<int>,int> &over_channel);
    double changing_obj();
    void find_violate();
    void create_pathGraph(vector<vector<vector<double>>> &path_graph_X, vector<vector<vector<double>>> &path_graph_Y, vector<vector<double>> &X_value, vector<vector<double>> & Y_value,bool print);

    MyAlgo5(string filename, int request_time_limit, int node_time_limit, double swap_prob, double entangle_alpha, int path_num);
    MyAlgo5(string filename, int request_time_limit, int node_time_limit, double swap_prob, double entangle_alpha, double epsilon , int path_num);
    ~MyAlgo5();
};

#endif
