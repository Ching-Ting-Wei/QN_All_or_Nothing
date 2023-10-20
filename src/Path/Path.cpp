#include "Path.h"

Path::Path(vector<Node*> nodes, vector<Channel*> channels)
    :nodes(nodes), channels(channels), swap_status(PATH_UNSWAP), send_status(PATH_UNSEND){
        prob = 1;
        for(auto channel:channels){
            prob *= channel->get_entangle_prob();
        } 
        for(int i = 1; i < (int)nodes.size()-1; i++){ // source and destination no swap
            prob *= nodes[i]->get_swap_prob();
        }

        //check if the path is connected
        Node *node1_ptr = nullptr, *node2_ptr = nodes[0];
        for(int i = 1; i < (int)nodes.size(); i++){
            node1_ptr = node2_ptr;
            node2_ptr = nodes[i];
            if(channels[i-1]->get_node1_ptr() == node1_ptr && channels[i-1]->get_node2_ptr() == node2_ptr)    continue;
            if(channels[i-1]->get_node2_ptr() == node1_ptr && channels[i-1]->get_node1_ptr() == node2_ptr)    continue;
            cerr << "error:\t the node and channel in path are not connected!" << endl;
            exit(1);
        }
        if(DEBUG)cerr << "new Path" << endl;
}

// Path::Path(Path & oldpath){
//     if(DEBUG)cerr << "new Path[COPY]" << endl;
// };

// Path::Path(vector<Node*> nodes)
//     :nodes(nodes), entangle_succ(true), swap_succ(false){
//         cerr<<"new entangled path"<<endl;
// }

Path::~Path(void){
    if(DEBUG)cerr<<"delete Path"<<endl;
    //release();
}

int Path::get_len(){
    return channels.size();
}

double Path::get_prob(){
    return prob;
}

vector<Node*> Path::get_nodes(){
    return nodes;
}

vector<Channel*> Path::get_channels(){
    return channels;
}

void Path::entangle(){                                  //question:why checking is_entangled do not put in this?
    if(DEBUG)cerr << "entangle path :" << endl;
    // this->print();
    if(channels.size() == 0 || nodes.size() <= 1){
        cerr << "error:\tWTF" << endl;
        exit(1);
    }
    bool entangle_succ = true;
    for(auto &channel : channels){
        if(DEBUG)cerr<<"entangle channel : ";
        entangle_succ &= channel->entangle();
    }

}

bool Path::get_entangle_succ(){
    bool entangle_succ = true;
    for(auto &channel : channels){
        if(DEBUG)cerr<<"entangle channel : ";
        entangle_succ &= channel->is_entangled();
    }
    return entangle_succ;
}

bool Path::swap(){
    if(!get_entangle_succ()){
        cerr<<"error:\tcan't swap on a path with fail entanglement!"<<endl;
        exit(1);
    }
    if(swap_status != PATH_UNSWAP){
        cerr<<"error:\tswap twice on the same path"<<endl;
        exit(1);
    }
    bool swap_succ = true;
    for(int i = 1; i < (int)nodes.size() - 1; i++){
        Node* node = nodes[i];
        swap_succ &= node->swap();
    }
    if(swap_succ){
        swap_status = PATH_SWAP_SUCC;
    }else{
        swap_status = PATH_SWAP_FAIL;
    }
    return swap_succ;
}

bool Path::get_swap_succ(){
    if(swap_status == PATH_UNSWAP){
        if(DEBUG)cerr<<"warning\tthe path is unswaped!"<<endl;
    }
    return (swap_status == PATH_SWAP_SUCC);
}

bool Path::send_data(){
    if(!get_swap_succ()){
        cerr<<"error:\tcan't send data on a unswaped path!"<<endl;
        exit(1);
    }
    if(send_status != PATH_UNSEND){
        cerr<<"error:\tthe path have already send the data!"<<endl;
        exit(1);
    }
    if(channels.size() == 0 || nodes.size() <= 1){
        cerr << "error:\tWTF(send data)" << endl;
        exit(1);
    }
    bool send_succ = true;
    for(auto &channel:channels){
        if(DEBUG)cerr<<"send data on channel : ";
        send_succ &= channel->send();
    }
    if(send_succ){
        send_status = PATH_SEND_SUCC;
    }else{
        send_status = PATH_SEND_FAIL;
    }
    return send_succ;
}

bool Path::send_data_succ(){
    if(!get_swap_succ()){
        cerr<<"error:\tcan't send data on a unswaped path!"<<endl;
        exit(1);
    }
    if(send_status == PATH_UNSEND){
        cerr<<"error:\tthe path haven't send the data yet!"<<endl;
        exit(1);
    }
    return (swap_status == PATH_SEND_SUCC);
}

void Path::release(){
    if(DEBUG)cerr<<"release()"<<endl;
    for(auto &channel:channels){
        channel->release();
    }
    send_status = PATH_UNSEND;
    swap_status = PATH_UNSWAP;
}
void Path::refresh(){
    for(auto &channel:channels){
        channel->refresh();
    }
    send_status = PATH_UNSEND;
    swap_status = PATH_UNSWAP;
}
void Path::print(){
    if(!DEBUG)  return;
    cerr<<"this is a Path"<<endl;
    for(auto n:nodes){
        cerr<<n->get_id()<<" ";
    }
    cerr << endl;
    cerr<<"\tentangle_succ = "<<get_entangle_succ()<<endl;
    cerr<<"\tswap_succ = "<<swap_status<<endl;
    cerr<<"\t the node in path is:";
    // for(auto n:nodes){
    //     cerr<<n->get_id()<<" ";
    // }
    // cerr << endl;
    cerr << "\t the channel in path is:";
    for(auto c : channels){
        cerr << c <<" ";
    }
    cerr << endl;
}

double Path::fidelity(){
    double fidelity = 1;
    /*                                          question:we just kill it ><
    for(auto &channel:channels){
        fidelity *= channel->get_fidelity();
    }
    */
    return fidelity;
}

bool Path::operator==(const Path &r)const{
    if(nodes.size() != r.nodes.size()) return false;
    for(int i=0;i<(int)nodes.size();i++){
        if(nodes[i] != r.nodes[i]) return false;
    }
    return true;
}

bool Path::operator<(const Path &r)const{
    if(nodes.size() != r.nodes.size()) return nodes.size() < r.nodes.size();
    for(int i = 0; i < (int)nodes.size(); i++){
        if(nodes[i] != r.nodes[i])  return nodes[i]->get_id() < r.nodes[i]->get_id();
    }
    return nodes.back() == r.nodes.back();
}