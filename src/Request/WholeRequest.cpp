#include "WholeRequest.h"

WholeRequest::WholeRequest(int source, int destination, int time_limit, vector<int> trusted_node_path)
    :Request(source, destination, time_limit), current_temporary(0), encode_to_5_qubits(false),
     encode_cnt(0), finished_qubits(0), success_qubits(0), path_length(0), fidelity(1), trusted_node_path(trusted_node_path){
    if(DEBUG)cerr<<"new WholeRequest"<<endl;
}

double calculate_fidelity(vector<double> paths_fidelity) {
    assert(!paths_fidelity.empty());
    sort(paths_fidelity.rbegin(), paths_fidelity.rend());
    vector<double> fidelities;
    while(fidelities.size() < 5LL) {
        for(int i = 0; i < (int)paths_fidelity.size() && fidelities.size() < 5LL; i++) {
            fidelities.push_back(paths_fidelity[i]);
        }
    }

    double product = 1;
    for(double fidelity : fidelities) {
        product *= fidelity;
    }

    double sum = product;
    for(double fidelity : fidelities) {
        sum += (product / fidelity) * (1 - fidelity);
    }

    return sum;
}

double WholeRequest::calculate_subfidelity() {
    if(subrequests.empty()) return 1;
    if(!is_encoding()) return subrequests[0]->get_paths()[0]->fidelity();

    vector<double> fidelities;
    for(SubRequest* subrequest : subrequests) {
        if(subrequest->get_paths().size() > 1) {
            cerr << "error subrequest has more than 1 path" << endl;
        }

        assert(!subrequest->get_paths().empty());
        fidelities.push_back(subrequest->get_paths()[0]->fidelity());
    }

    return calculate_fidelity(fidelities);
}

int WholeRequest::get_current_temporary() {
    return current_temporary;
}
int WholeRequest::get_trust_node_path_length(){
    return trusted_node_path.size();
}
map<int, int> WholeRequest::get_use_path_count(){
    return use_path_count;
}
int WholeRequest::get_encode_cnt(){
    return encode_cnt;
}

bool WholeRequest::is_encoding() {
    return encode_to_5_qubits;
}
bool WholeRequest::is_finished() {
    return current_temporary == (int)trusted_node_path.size() - 1;
}

bool WholeRequest::is_success() {
    return status == REQUEST_SUCC;
}

void WholeRequest::set_encode(bool flag) {
    encode_to_5_qubits = flag;
}

void WholeRequest::temporary_forward() {
    current_temporary++;
    fidelity *= calculate_subfidelity();
    if(is_finished()) {
        if(status == REQUEST_UNFINISHED) {
            status = REQUEST_SUCC;
        }
    }
    Node *dst = subrequests.back()->get_paths()[0]->get_nodes().back();
    for(SubRequest* subrequest : subrequests) {
        delete subrequest;
        (*dst)--;
    }

    if(is_finished()) {
        (*dst)++;
    }

    if(encode_to_5_qubits) {
        (*dst) += 4;
    }

    encode_to_5_qubits = false;
    finished_qubits = 0;
    success_qubits = 0;
    subrequests.clear();
}

void WholeRequest::try_forward() {
    for(int i=0;i<(int)subrequests.size();i++) {
        SubRequest* subrequest = subrequests[i];
        if(subrequest->is_finished()) {
            path_count.insert(i);
            if(is_encoding()) path_length += subrequest->get_paths()[0]->get_len() / (double)5;
            else path_length += subrequest->get_paths()[0]->get_len();
            finished_qubits++;
            Node *src = subrequest->get_paths()[0]->get_nodes()[0];
            (*src)++;
            if(subrequest->is_success()) {
                success_qubits++;
            }
        }
    
    }

    if(is_encoding()) {                           //divide into five qubits
        if(finished_qubits >= 5) {
            if(finished_qubits > 5) {
                cerr << "error: finished_qubits > 5" << endl;
            }
            if(status != REQUEST_FAIL && success_qubits >= 4) {
                status = REQUEST_UNFINISHED;
            } else {
                status = REQUEST_FAIL;
            }

            encode_cnt++;
            use_path_count[path_count.size()]++;
            path_count.clear();
            temporary_forward();
        } else {
            int need = (5 - finished_qubits);
            int remove_cnt = (int)subrequests.size() - need;
            for(int i = 0; i < remove_cnt; i++) {
                Node *dst = subrequests.back()->get_paths()[0]->get_nodes().back();
                delete subrequests.back();
                subrequests.pop_back();
                (*dst)--;
            }
        }
    } else {
        if(finished_qubits >= 1) {
            if(finished_qubits > 1) {
                cerr << "error: finished_qubits > 1, but undivide" << endl;
            }

            if(status != REQUEST_FAIL && success_qubits == 1) {
                status = REQUEST_UNFINISHED;
            } else {
                status = REQUEST_FAIL;
            }
            temporary_forward();
        }
    }
}

void WholeRequest::clear() {
    for(SubRequest* subrequest : subrequests) {
        delete subrequest;
    }
}