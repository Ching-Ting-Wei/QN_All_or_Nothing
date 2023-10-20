#include "Greedy.h"

Greedy::Greedy(string filename, int request_time_limit, int node_time_limit, double swap_prob, double entangle_alpha , bool limit_r_or_not)
    :AlgorithmBase(filename, "Greedy", request_time_limit, node_time_limit, swap_prob, entangle_alpha , limit_r_or_not){
    if(DEBUG) cerr<<"new Greedy"<<endl;
}

void Greedy::path_assignment(){
    // base_test_active();
    bool flag = true;
    while(flag){
        flag = false;
        for(int reqno = 0; reqno < (int)requests.size(); reqno++){
            Request &request = requests[reqno];
            if(get_limit_r_status()){
                if(request.get_paths().size() >= request.get_send_limit()){
                    continue;
                }
            }
            vector<int> path = BFS(request.get_source(), request.get_destination());
            if(path.size() == 0 ) continue;
            assign_resource(path, reqno);
            flag = true;
        }
    }
}

void Greedy::entangle(){
    AlgorithmBase::base_entangle();
}
void Greedy::swap(){
    AlgorithmBase::base_swap();
}
void Greedy::send(){
    AlgorithmBase::base_send();
}
void Greedy::next_time_slot(){
    AlgorithmBase::base_next_time_slot();
}