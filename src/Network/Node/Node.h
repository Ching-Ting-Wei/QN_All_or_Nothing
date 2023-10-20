#ifndef __NODE_H
#define __NODE_H

#include <iostream>
#include <random>
#include <ctime>
#include <cmath>
#include "../../config.h"
using namespace std;

class Node{
    int id;
    double swap_prob;
    int memory_cnt, time_limit;
    int remain;                 //remain qubit
    pair <double, double> pos;   //唯一決定一個 node
    void release();                              //clean assign status
public:
    Node(int id, int memory_cnt, int time_limit, double pos_x, double pos_y, double swap_prob);
    // Node(Node & old_node);
    //test
    ~Node();
    Node():swap_prob(0){}
    int get_id();
    int get_memory_cnt();
    double get_swap_prob();

    bool swap();                 //release one remain memory
    double distance(const Node &right)const;     //return the distance of two node
    bool is_assignable()const;                   //return if we can use this node to build an entangle or not
    int get_remain()const;                       //return the number of qubit is not used 
    void print()const;
    void revise(int mem);
    bool operator==(const Node &right)const;
    bool operator!=(const Node &right)const;
    bool operator<(const Node &right)const;
    bool operator<=(const Node &right)const;
    bool operator>(const Node &right)const;
	bool operator>=(const Node &right)const;
    const Node operator--(int);                  //delete one remain memory
    const Node operator++(int);
    const Node operator+=(int);
    const Node operator-=(int);
};

#endif