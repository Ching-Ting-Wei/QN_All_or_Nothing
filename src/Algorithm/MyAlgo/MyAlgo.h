#ifndef __MYALGO_H
#define __MYALGO_H

#include <iostream>
#include <algorithm>
#include <omp.h>
#include <queue>
#include <limits>
#include <string>
#include <cmath>
#include "../AlgorithmBase/AlgorithmBase.h"
#include "../../Request/WholeRequest.h"
#include "../../Network/Graph/Graph.h"
#include "../../config.h"

using namespace std;

class MyAlgo:public AlgorithmBase {

private:
    map<pair<int,int>, vector<double>> X;
    vector<map<pair<int,int>, double>> Y;
    vector<double> alpha;
    map<vector<int>, double> x_i_p;
    map<pair<int,int>, double> beta;
    vector<vector<vector<int>>> all_source_target_path;
    int change_edge_num = 0; 
    double epsilon = 0.2;
    double delta;
    double M;
    vector<double> tau;                
public: 
    map<int, int> num_of_path_count;
    map<int, int> path_length_encode;
    map<int, int> path_length_cnt;
    vector<int> separation_oracle(int req_no, double &U);
    vector<int> Dijkstra(int src, int dst, int req_no);
    vector<int> Dijkstra_by_prob(int src, int dst);
    void path_assignment();
    void calculate();
    void entangle();
    void swap();
    void send();
    void dfs(int s, int t, vector<vector<int>> &ans, vector<int> &path, vector<bool> &visited);    
    vector<vector<int>> allPathsSourceTarget(int src, int dst);
    void next_time_slot();
    vector<map<vector<int>, int>> Greedy_rounding();
    void find_bottleneck(vector<int>, int req_no);
    void initialize();
    void check_enough(vector<map<vector<int>, int>> &path);
    void readd(vector<map<vector<int>, int>> &path,vector<int> &over_memory,map<vector<int>,int> &over_channel);
    double changing_obj();
    void find_violate();
    vector<map<vector<int>, int>> rounding();
    MyAlgo(string filename, int request_time_limit, int node_time_limit, double swap_prob, double entangle_alpha);
    ~MyAlgo();
};

#endif
