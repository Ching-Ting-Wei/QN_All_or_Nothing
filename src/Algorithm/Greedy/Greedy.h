#ifndef __GREEDY
#define __GREEDY
#include <iostream>
#include <algorithm>
#include <queue>
#include "../AlgorithmBase/AlgorithmBase.h"
#include "../../Network/Graph/Graph.h"
#include "../../config.h"
using namespace std;

class Greedy:public AlgorithmBase{
public:
    Greedy(string filename, int request_time_limit, int node_time_limit, double swap_prob, double entangle_alpha , bool limit_r_or_not);
    void path_assignment();
    void entangle();
    void swap();
    void send();
    void next_time_slot();
};

#endif