#include"Request.h"

Request::Request(int source, int destination, const int& time_limit):source(source),
    destination(destination), time_limit(time_limit), status(REQUEST_UNFINISHED), send_path_length(0){
    if(DEBUG)cerr<<"new Request"<<endl;
}

Request::Request(int source, int destination, const int& time_limit, int send_limit):source(source),
    destination(destination), time_limit(time_limit), send_limit(send_limit), status(REQUEST_UNFINISHED), send_path_length(0){
    if(DEBUG)cerr<<"new Request"<<endl;
}

Request::~Request(void){
    if(DEBUG)cerr<<"delete Request"<<endl;
    for(int i=0;i<(int)paths.size();i++){
        paths[i]->release();
        delete paths[i];
        paths[i] = nullptr;
    }
    paths.clear();
}

void Request::set_path(int path_id, Path *p){
    if(path_id >= (int)paths.size()){
        cerr<<"can't set this path!"<<endl;
        exit(1);
    }
    paths[path_id] = p;
}

int Request::get_time_limit(){
    return time_limit;
}

int Request::get_send_limit(){
    return send_limit;
};

int Request::get_path_num(){
    return path_num;
};

int Request::get_before_ent_path_num(){
    return before_ent_path_num;
};

double Request::get_before_ent_total_prob(){
    return before_ent_total_prob;
}
// int Request::get_waiting_time(){
//     return waiting_time;
// }
int Request::get_source(){ 
    return source;
}
int Request::get_destination(){
    return destination;
}
int Request::get_send_path_length(){
    if(status == REQUEST_UNFINISHED){
        cerr<<"error:\tthe request is unfinished!"<<endl;
        exit(1);
    }
    return send_path_length;
}

double Request::get_total_prob(){
    return total_prob;
}

double Request::get_fidelity(){
    if(status == REQUEST_UNFINISHED){
        cerr<<"error:\tthe request is unfinished!alkdjf"<<endl;
        exit(1);
    }
    return fidelity;
}
vector<Path *> Request::get_paths(){
    return paths;
}

int Request::get_throughput(){
    return this->throughput;
}

void Request::clear_paths(){
    for(int i=0;i < (int)paths.size();i++){
        if(paths[i] != nullptr){
            delete paths[i];
            paths[i] = nullptr;
        }
    }
    paths.clear();
}

void Request::refresh_paths(){
    for(Path *path : paths){
        path->refresh();
    }
    status = REQUEST_UNFINISHED;
}

void Request::entangle(){
    before_ent_path_prob_vt.clear();
    for(auto &path:paths){
        if(path == nullptr)continue;
        before_ent_total_prob += path->get_prob();
        before_ent_path_num++;
        before_ent_path_prob_vt.push_back(path->get_prob());
        path->entangle(); 
    }
}

void Request::swap(){
    cout<< "request limit: " << send_limit << endl;
    cout<<"swapping path number: " << paths.size() << endl;
    for(auto &path:paths){
        if(path == nullptr)continue;
        total_prob += path->get_prob();
        path_num++;
        if(path->get_entangle_succ()) {  
            if( path->swap()){
                throughput++;
            }
        }
        
    }
}

void Request::send(){
    int pid = -1, mx = 0;
    for(int i=0;i<(int)paths.size();i++){
        if(paths[i] == nullptr)continue;
        if(!paths[i]->get_swap_succ())continue;
        if(paths[i]->fidelity() > mx){
            mx = paths[i]->fidelity();
            pid = i;
        }
    }
    if(pid == -1){
        return;
    }

    send_path_length = paths[pid]->get_len();
    fidelity = paths[pid]->fidelity();
    if(paths[pid]->send_data()){
        status = REQUEST_SUCC;
    }else{
        status = REQUEST_FAIL;
    }
}
bool Request::is_finished(){
    return status != REQUEST_UNFINISHED;
}
bool Request::is_success(){
    if(status == REQUEST_UNFINISHED){
        cerr<<"the request is unfinshed!"<<endl;
        exit(1);
    }
    return status == REQUEST_SUCC;
}
void Request::next_timeslot(){
    for(auto path_ptr:paths){
        if(path_ptr != nullptr){
            path_ptr->release();
            delete path_ptr;
        }
    }
    paths.clear();
    // waiting_time++;
    // if(throughput > 0){
    //     time_limit--;
    // }
    // return (time_limit == 0) && (throughput > 0);
}

void Request::add_cur(int num){
    cur_send += num;
}

int Request::get_cur_send(){
    return cur_send;
}

vector<double> Request::get_before_ent_path_prob_vt(){
    return before_ent_path_prob_vt;
}

void Request::operator+=(Path *path){
    paths.emplace_back(path);
}

void Request::delete_path(){
    paths.pop_back();
}



