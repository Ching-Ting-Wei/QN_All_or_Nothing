#include "AlgorithmBase.h"



AlgorithmBase::AlgorithmBase(string filename, string algorithm_name, int request_time_limit, int node_time_limit, double swap_prob, double entangle_alpha , bool limit_r_or_not)
    :timeslot(0), waiting_time(0), time_limit(request_time_limit), swap_prob(swap_prob), limit_r_or_not(limit_r_or_not) , graph(Graph(filename, node_time_limit, swap_prob, entangle_alpha) ){
    if(limit_r_or_not==true){
        this->algorithm_name=algorithm_name;
    }
    else{
       this->algorithm_name=algorithm_name+"_Nonlimit";
    }
}

AlgorithmBase::~AlgorithmBase(){
    if(DEBUG)cerr<<"delete AlgorithmBase"<<endl;
}

double AlgorithmBase::get_swap_prob(){
    return swap_prob;
}


bool AlgorithmBase::get_limit_r_status(){
    return limit_r_or_not;
}

void AlgorithmBase::base_next_time_slot(){
    int total_path_num = 0;
    int before_ent_path_num = 0;
    double total_success_prob = 0;
    double before_ent_total_success_prob = 0;
    double min_req_success_ratio = 100;
    double max_req_success_ratio = 0;
    double total_req_success_ratio = 0;
    double max_over_ratio = 0;
    int total_earn = 0;
    int drop_req_no = 0;
    graph.refresh();
    graph.release();
    for(auto &request: requests){
        request.next_timeslot();
    }

    //好強
    vector<int> finished_reqno;
    res_vt.clear();
    for(int reqno = 0; reqno < (int)requests.size(); reqno++) {
        double req_success_ratio;
        total_path_num += requests[reqno].get_path_num();
        before_ent_path_num += requests[reqno].get_before_ent_path_num();
        total_success_prob += requests[reqno].get_total_prob();
        before_ent_total_success_prob += requests[reqno].get_before_ent_total_prob();

        for(auto it:requests[reqno].get_before_ent_path_prob_vt()){
            res_vt.push_back(it);
        }
        if(requests[reqno].get_path_num() >= requests[reqno].get_send_demand()){
            total_earn += requests[reqno].get_value();
        }
        else{
            drop_req_no ++;
        }
        if(requests[reqno].get_throughput() == 0){
            req_success_ratio = 0.0;
        }else{
            req_success_ratio = (double)requests[reqno].get_throughput() / (double)requests[reqno].get_send_demand();
            if(max_over_ratio <= req_success_ratio){
                max_over_ratio = req_success_ratio;
            }
            if(req_success_ratio >= 1){
                req_success_ratio = 1;
            }
        }

        if(min_req_success_ratio > req_success_ratio){
            min_req_success_ratio = req_success_ratio;
        }
        if(max_req_success_ratio < req_success_ratio){
            max_req_success_ratio = req_success_ratio;
        }
        total_req_success_ratio += (1 - req_success_ratio);
        if(!requests[reqno].is_finished()) {
            continue;
        }
        res["finished_throughputs"]++;
        res["path_length"] += requests[reqno].get_send_path_length();
        res["fidelity"] += requests[reqno].get_fidelity();
        
        finished_reqno.push_back(reqno);
        if(requests[reqno].get_throughput() != 0){
            cout << "reqno: " << reqno << " " <<  requests[reqno].get_throughput() << endl;
            if(requests[reqno].get_throughput() >= requests[reqno].get_send_demand()){
                res["throughputs"]+= requests[reqno].get_send_demand();
            }else{
                res["throughputs"]+= requests[reqno].get_throughput();
            }  
        }
    }
    res["total_earn"] = total_earn;
    res["drop_req_no"] = drop_req_no;
    res["path_success_avg"] = total_success_prob / total_path_num;
    res["path_success_avg_before_ent"] = before_ent_total_success_prob / before_ent_path_num;
    res["S_D_complete_ratio_difference"] = max_req_success_ratio - min_req_success_ratio;
    res["new_success_ratio"] = total_req_success_ratio / requests.size();
    res["max_over_ratio"] = max_over_ratio;
    res["dis_avg"] = graph.get_dis_avg();
    res["prob_avg"] = graph.get_prob_avg();
    reverse(finished_reqno.begin(), finished_reqno.end());
    for(int reqno : finished_reqno) {
        requests.erase(requests.begin() + reqno);
    }
   
}

void AlgorithmBase::base_entangle(){
    for(int i=0;i<(int)requests.size();i++){
        requests[i].entangle();
    }
}

void AlgorithmBase::base_swap(){
    for(int i=0;i<(int)requests.size();i++){
        requests[i].swap();
    }
}
void AlgorithmBase::base_send(){
    for(int i=0;i<(int)requests.size();i++){
        requests[i].send();
   }
}

void AlgorithmBase::base_test_active(){
    if(requests.size() != 0) res["active_timeslot"]++;
    res["waiting_time"] += requests.size();
}

void AlgorithmBase::check_resource(){
    double use_memory_sum = 0, memory_sum = 0, use_channel_sum = 0, channel_sum = 0;
    for(int i = 0; i < (int)graph.get_size(); i++){
        Node *node_ptr = graph.Node_id2ptr(i);
        use_memory_sum += node_ptr->get_memory_cnt() - node_ptr->get_remain();
        memory_sum += node_ptr->get_memory_cnt();
    }
    res["use_memory"] = max(res["use_memory"], use_memory_sum);
    res["total_memory"] = memory_sum;
    res["memory_cnt_avg"] = memory_sum/(int)graph.get_size();
    for(int i = 0; i < (int)graph.get_size(); i++){
        for(int j = i+1; j < (int)graph.get_size(); j++){
            use_channel_sum += graph.get_used_channel_size(i, j);
            channel_sum += graph.get_channel_size(i, j);
        }
    }
    res["use_channel"] = max(res["use_channel"], use_channel_sum);
    res["total_channel"] = channel_sum;
}

void AlgorithmBase::run(){
    timeval start;
    timeval end;
    gettimeofday(&start, NULL);
    path_assignment();
    check_resource();
    entangle();
    swap();
    send();
    cout<<algorithm_name<<" run end"<<endl;
    next_time_slot();
    gettimeofday(&end, NULL);
    res["runtime"] += (1000000LL * (end.tv_sec - start.tv_sec) + end.tv_usec - start.tv_usec) / 1000000.0;
}

void AlgorithmBase::add_new_request(Request new_request){
    int src = new_request.get_source();
    requests.push_back(new_request);
    // if(graph.Node_id2ptr(src)->get_remain() >= 1) {
    //     requests.push_back(new_request);
    //     (*graph.Node_id2ptr(src))--;
    // }
}

int AlgorithmBase::find_width(vector<int> path){
    if(path.size() < 2){
        return 0x3f3f3f3f;
    }
    if(path.size() == 2){
        return graph.remain_resource_cnt(path[0], path[1], false, false);
    }
    int path_size = (int)path.size();
    int width = graph.remain_resource_cnt(path[0], path[1], false);
    for(int i = 2; i < path_size - 1; i++){
        width = min(width, graph.remain_resource_cnt(path[i-1], path[i]));
    }
    width = min(width, graph.remain_resource_cnt(path[path_size-2], path[path_size-1], true, false));
    return width;
}

void AlgorithmBase::assign_resource(vector<int> path, int reqno){
    
    if(DEBUG) cerr<< "---------AlgorithmBase::assign_resource----------" << endl;
    int width = find_width(path);
    if(DEBUG) {
        cerr<<"find a path to build! with width:    "<<width<<endl;
        for(auto ele:path){
            cerr<<ele<<" ";
        }
        cerr << endl;
    }
    if(width == 0x3f3f3f3f){
        cerr<<"error:\twidth = INF"<<endl;
        exit(1);
    }
    while(width-- > 0){
       
        requests[reqno] +=  graph.build_path(path);
    }
    if(DEBUG) cerr<< "---------AlgorithmBase::assign_resource----------end" << endl;
}

void AlgorithmBase::assign_resource(vector<int> path, int width, int reqno){
    
    if(DEBUG) cerr<< "---------AlgorithmBase::assign_resource----------" << endl;
    if(DEBUG) {
        cerr<<"find a path to build! with width:    "<<width<<endl;
        for(auto ele:path){
            cerr<<ele<<" ";
        }
        cerr << endl;
    }
    if(width == 0x3f3f3f3f){
        cerr<<"error:\twidth = INF"<<endl;
        exit(1);
    }

    while(width-- > 0){
       
        requests[reqno] +=  graph.build_path(path);
    }
    if(DEBUG) cerr<< "---------AlgorithmBase::assign_resource----------end" << endl;
}

vector<int> AlgorithmBase::BFS(int source, int destination){
    // if(source > graph.get_size() || destination > graph.get_size()){
    //     cerr<<"error:\tnode id is out of range"<<endl;
    //     exit(1);
    // }
    vector<bool> vis(graph.get_size()+5, false);
    vector<int> parent(graph.get_size()+5, -1);
    // for(int i = 0; i < graph.get_size(); i++){
    //     vis[i] = false;
    //     parent[i] = -1;
    // }

    queue<int> q;
    q.push(source);
    vis[source] = true;
    while (!q.empty()){
        int now = q.front();
        q.pop();
        vector<int> neighbors = graph.get_neighbors_id(now);
        for(auto neighbor: neighbors){
            if(vis[neighbor]) continue;
            bool is1_repeater = true, is2_repeater = true;
            if(now == source) is1_repeater = false;
            if(neighbor == destination) is2_repeater = false;
            if(graph.remain_resource_cnt(now, neighbor, is1_repeater, is2_repeater) > 0){
                parent[neighbor] = now;
                q.push(neighbor);
                vis[neighbor] = true;
            }
        }
    }
    vector<int> path_nodes;
    path_nodes.clear();
    if(parent[destination] == -1) return path_nodes;

    int now = destination;
    while(now != -1){
        path_nodes.push_back(now); 
        now = parent[now];
    }
    reverse(path_nodes.begin(), path_nodes.end());
    return path_nodes;
}

vector<Request> AlgorithmBase::get_requests(){
    return requests;
};

int AlgorithmBase::total_throughput(){
    return res["throughputs"];
}

Path* AlgorithmBase::find_swap_path(vector<int> path_nodes, map<pair<int, int>, vector<Channel*>> &remain_channels){
    //the the path for swaping
    //cout << "AlgorithmBase::find_swap_path" << endl;
    
    
    //find the channels of path for swaping 
    vector<Channel*> path_channels;
    int u, v;
    for(int i=1;i<(int)path_nodes.size();i++){
        u = path_nodes[i-1];
        v = path_nodes[i];
        if(u > v){
            int swap_tmp = u;
            u = v;
            v = swap_tmp;
        }
        if(remain_channels[make_pair(u, v)].empty()){
            cerr<<"error:\t remain channel is not sufficient"<<endl;
            exit(1);
        }
        path_channels.emplace_back(remain_channels[make_pair(u, v)].back());
        remain_channels[make_pair(u, v)].pop_back();
    }

    //find the nodes of path for swapin
    vector<Node *> nodes;
    if(DEBUG) cerr << "path to swap in P4: ";
    for(int i:path_nodes){
        nodes.emplace_back(graph.Node_id2ptr(i));
        if(DEBUG) cerr <<  i << ' ';
    }
    if(DEBUG) cerr << '\n';
    return new Path(nodes, path_channels);
}

string AlgorithmBase::get_name(){
    return algorithm_name;
}

double AlgorithmBase::get_res(string s){
    return res[s];
}

vector<double> AlgorithmBase::get_res_vt(){
    return res_vt;
};

