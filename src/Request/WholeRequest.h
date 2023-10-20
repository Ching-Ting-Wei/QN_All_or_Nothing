#ifndef __WHOLEREQUEST_H
#define __WHOLEREQUEST_H
#include    <iostream>
#include    <vector>
#include    <algorithm>
#include    <cassert>
#include    <set>
#include    <map>
#include    "../config.h"
#include    "Request.h"
using namespace std;
using SubRequest = Request;                 //genius


class WholeRequest:public Request {
    int current_temporary; // index of node in trusted_node_path
    bool encode_to_5_qubits;
    map<int, int> use_path_count;
    set<int> path_count;
    set<int> node_path_count;
    map<int, int> use_node_path_count;
    int encode_cnt;
public:
    int finished_qubits;
    int success_qubits;
    double path_length;
    double fidelity;
    vector<int> trusted_node_path;
    vector<SubRequest*> subrequests;
    WholeRequest(int source, int destination, int time_limit, vector<int> trusted_node_path);
    double calculate_subfidelity();
    int get_current_temporary();
    int get_trust_node_path_length();
    int get_encode_cnt();
    map<int, int> get_use_path_count();
    bool is_encoding();
    bool is_finished();
    bool is_success();
    void set_encode(bool status);
    void try_forward();
    void temporary_forward();
    void clear();
};

#endif