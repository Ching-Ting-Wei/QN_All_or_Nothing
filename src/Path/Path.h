#ifndef __PATH_H
#define __PATH_H

#include <vector>
#include "../Network/Node/Node.h"
#include "../Network/Channel/Channel.h"
#include "../config.h"

const int PATH_UNSWAP = 0;
const int PATH_SWAP_SUCC = 1;
const int PATH_SWAP_FAIL = -1;

const int PATH_UNSEND = 0;
const int PATH_SEND_SUCC = 1;
const int PATH_SEND_FAIL = -1;


class Path{
    vector<Node*> nodes;
    vector<Channel*> channels;
    int swap_status;
    int send_status;
    double prob;
public:
	Path(vector<Node*> nodes, vector<Channel*> channels);
    // Path(Path & oldpath);
    //Path(vector<Node*> nodes);                        zzzz
    ~Path(void);
    int get_len();
    vector<Node*> get_nodes();
    vector<Channel*> get_channels();
    void entangle();        //true->entangle success
    bool get_entangle_succ();
    bool swap();            //true->swapped success
    bool get_swap_succ();
    bool send_data();
    bool send_data_succ();
    void release();
    void refresh();
    void print();
    double fidelity();
    double get_prob();
    bool operator==(const Path &r)const;
    bool operator<(const Path &r)const;
};

#endif