#include "../AlgorithmBase/AlgorithmBase.h"
#include "../../Network/Graph/Graph.h"
#include "../../config.h"
#include "gurobi_c++.h"
using namespace std;


class REPS:public AlgorithmBase{
    vector<Channel *> assign_path;
    //vector<double> t_hat;
    //vector<map<pair<int, int>, int>> f_hat;
    vector<vector<vector<int>>> EPS_P;
    vector<vector<vector<int>>> ELS_P;
    void PFT_LP(vector<double> &t_plum, vector<map<pair<int, int>, double>> &f_plum);
    void EPS_LP(vector<vector<double>> &t_bar, vector<vector<map<pair<int, int>, double>>> &f_bar);
    void EPS();
    void ELS();
    tuple<vector<int>, double, bool> DFS(int req_no, map<pair<int, int>, double>&f_plum_i, bool is_assign_path = true);
public:
    REPS(string filename, int request_time_limit, int node_time_limit, double swap_prob, double entangle_alpha , bool limit_r_or_not);
    void path_assignment();
    void entangle();
    void swap();
    void send();
    void next_time_slot();
};