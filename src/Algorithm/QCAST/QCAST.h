#include    <iostream>
#include    <cmath>
#include    <algorithm>
#include    <stack>
#include    "../AlgorithmBase/AlgorithmBase.h"
#include    "../../Network/Graph/Graph.h"
#include "../../config.h"
using namespace std;

/* 
找 major path 的方法: Greedy-EDA(G-EDA)
Step1)  用 Extended Dijkstra's algorithm 根據權重 EXT 找出每個 S-D pair 的最佳 path
Step2)  選最大 EXT 的 path，分配資源後更新為 residual graph
重複 step 1, 2
直到無法再找到 path，或找到的 path 超過一定數量
設定 path length 的最大值，超過就直接 pop 掉
(先不用recovery path)
*/ 
class QCAST:public AlgorithmBase{
    map<pair<Request*,Path*>,vector<Path*>> recovery_paths;
    double EXT(vector<double> path, int w);
    map<pair<int, int>, double> combination;
    double C(int n, int m);
    // void assign_resource(vector<int> path, int reqno);                   //assign path's resource to requests[reqno]
    void p4();
public:
    QCAST(string filename, int request_time_limit, int node_time_limit, double swap_prob, double entangle_alpha , bool limit_r_or_not);
    ~QCAST();
    void path_assignment();
    void find_recovery_path(int R);
    void entangle();
    void swap();
    void next_time_slot();
    void send();
    double demoEXT(vector<double> path, int w);
};