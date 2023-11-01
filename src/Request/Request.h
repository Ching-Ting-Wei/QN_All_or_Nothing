#ifndef __REQUEST_H
#define __REQUEST_H

#include<iostream>
#include<vector>
#include<algorithm>
#include"../Network/Node/Node.h"
#include"../Path/Path.h"
#include "../config.h"
using namespace std;

const int REQUEST_UNFINISHED = 0;
const int REQUEST_SUCC = 1;
const int REQUEST_FAIL = -1;


class Request{
protected:
    int source, destination;
    int time_limit;
    int send_demand;
    int value;
    int status;
    int send_path_length;

    int cur_send = 0;
    int throughput = 0;
    double willness=0;             //bug
    double fidelity;
    double total_prob = 0;
    int path_num = 0;
    double before_ent_total_prob = 0;
    int before_ent_path_num = 0;
    vector<Path *> paths;                                       //休學
    vector<double>before_ent_path_prob_vt;
public:
    Request(int source, int destination, const int& time_limit);
    Request(int source, int destination, const int& time_limit, int send_demand, int value);
    ~Request(void);
    void set_path(int path_id, Path *p);                        //should delete old path before set new path
    int get_time_limit();
    int get_send_demand();
    int get_value();
    int get_path_num();
    int get_before_ent_path_num();
    double get_before_ent_total_prob();

    int get_source();
    int get_destination();
    int get_send_path_length();
    double get_total_prob();
    double get_fidelity();
    vector<Path *> get_paths();
    int get_throughput();
    void clear_paths();
    void refresh_paths();
    void add_one_throughput();
    void entangle();
    void swap();
    void send();
    bool is_finished();
    bool is_success();
    void next_timeslot();
    void add_cur(int num);
    int get_cur_send();
    vector<double> get_before_ent_path_prob_vt();
    void operator+=(Path *path);
    void delete_path();

    double get_willness();            //bug
    void set_value(int val);          //bug
    void set_willness(double will);   //bug
    //void print();
    // int get_waiting_time();
};

#endif