#include "MyAlgo4.h"

MyAlgo4::MyAlgo4(string filename, int request_time_limit, int node_time_limit, double swap_prob, double entangle_alpha)
    :AlgorithmBase(filename, "MyAlgo4", request_time_limit, node_time_limit, swap_prob, entangle_alpha , true /*let it be alway true*/){
    if(DEBUG) cerr<<"new MyAlgo4"<<endl;
}
MyAlgo4::MyAlgo4(string filename, int request_time_limit, int node_time_limit, double swap_prob, double entangle_alpha, double epsilon)
    :AlgorithmBase(filename, "MyAlgo4_" + to_string(epsilon) , request_time_limit, node_time_limit, swap_prob, entangle_alpha , true /*let it be alway true*/), epsilon(epsilon){
    if(DEBUG) cerr<<"new MyAlgo4"<<endl;
}


MyAlgo4::~MyAlgo4(){
    if(DEBUG) cerr << "delete MyAlgo4" << endl;
}

void MyAlgo4::entangle(){
    AlgorithmBase::base_entangle();
}

void MyAlgo4::swap(){
     AlgorithmBase::base_swap();
}

void MyAlgo4::send(){
     AlgorithmBase::base_send();
}

void MyAlgo4::next_time_slot(){
     AlgorithmBase::base_next_time_slot();
}

double MyAlgo4::X(int u, int v, vector<double>&alpha , map<pair<int,int>,double>&beta ){
	double weight = alpha[u] + alpha[v] + beta[{u, v}];
	return weight;
}

void MyAlgo4::initialize(){
    M = graph.get_size() + graph.get_num_of_edge() + requests.size();              //M=V+E+I
    // delta = (1 + epsilon) * pow(((1 + epsilon) * M), (-1 / epsilon));         
    delta = (1 + epsilon) * (1.0 / pow((1 + epsilon) * M, 1.0 / epsilon));

    for(int i = 0; i < graph.get_size(); i++){
        alpha.emplace_back(delta / graph.Node_id2ptr(i)->get_memory_cnt());       //alpha=delta/Mu

        vector<int> temp = graph.get_neighbors_id(i);                             //beta=delta/Cuv
        for(auto it: temp){
            // if(i < it){
            //     beta[make_pair(i,it)] = delta / (graph.get_channel_size(i, it));
            // }
            // else{
            //     beta[make_pair(it,i)] = delta / (graph.get_channel_size(i, it));
            // }
            beta[make_pair(i,it)] = delta / (graph.get_channel_size(i, it));
            beta[make_pair(it,i)] = delta / (graph.get_channel_size(i, it));

        }
    }

    Y.resize(requests.size() + 5);
    for(int i = 0; i < graph.get_size(); i++){
        vector<int> temp = graph.get_neighbors_id(i);                             
        for(auto it: temp){

            for(unsigned  j = 0;j < requests.size(); j++){                       //Y=-ln(edge)-ln(repeater_1^(1/2))-ln(repeater_2^(1/2))
                int src = requests[j].get_source();
                int des = requests[j].get_destination();
                double ent_p = exp(graph.Node_id2ptr(i)->distance(*graph.Node_id2ptr(it))*(-graph.get_entangle_alpha()));

                if( make_pair(i, it) == make_pair(src, des) || make_pair(it, i) == make_pair(src, des) ) {
                    Y[j][{i, it}] = -log(ent_p);
                    Y[j][{it, i}] = -log(ent_p);
                }
                else if( i == src || i == des ){
                    Y[j][{i, it}] = -log(ent_p) - log(graph.Node_id2ptr(it)->get_swap_prob()) / 2;
                    Y[j][{it, i}] = -log(ent_p) - log(graph.Node_id2ptr(it)->get_swap_prob()) / 2;
                }
                else if( it == src || it == des){
                    Y[j][{i, it}] = -log(ent_p) - log(graph.Node_id2ptr(i)->get_swap_prob()) / 2;
                    Y[j][{it, i}] = -log(ent_p) - log(graph.Node_id2ptr(i)->get_swap_prob()) / 2;
                }
                else{
                    Y[j][{i, it}] = -log(ent_p) - log(graph.Node_id2ptr(i)->get_swap_prob() * graph.Node_id2ptr(it)->get_swap_prob()) / 2;
                    Y[j][{it, i}] = -log(ent_p) - log(graph.Node_id2ptr(i)->get_swap_prob() * graph.Node_id2ptr(it)->get_swap_prob()) / 2;
                }
            }
        }
    }
}

vector<int> MyAlgo4::Dijkstra(int src, int dst, int req_no){ 
    const double INF = numeric_limits<double>::infinity();
    //const double INF = 100;
    int n = graph.get_size();
    vector<double> distance(n, INF);
    vector<int> parent(n, -1);
    vector<bool> used(n, false);
    priority_queue<pair<double, int>, vector<pair<double, int>>, greater<pair<double, int>>> pq;

    distance[dst] = 0;
    pq.push({0, dst});
    while(!pq.empty()) {
        int cur_node = pq.top().second;
        pq.pop();
        if(used[cur_node]) continue;
        used[cur_node] = true;
        for(int neighbor : graph.get_neighbors_id(cur_node)) {
            if(distance[cur_node] + X(cur_node, neighbor, alpha , beta ) < distance[neighbor]) {
                distance[neighbor] = distance[cur_node] + X(cur_node, neighbor, alpha , beta );
                parent[neighbor] = cur_node;
                pq.push({distance[neighbor], neighbor});
            }
        }
    }

    if(distance[src] >= INF) return{};
    
    return parent;
}       
   
void MyAlgo4::combination(int index, int req_no, vector<int>&cur_comb){
    if(cur_comb.size()==3/*this is k*/){
        MyAlgo4::comb[req_no].push_back(cur_comb);
        return;
    }
    //cout<<" size():"<<all_source_target_path[req_no].size()<<endl;
    for(unsigned int i=index; i<MyAlgo4::all_source_target_path[req_no].size(); i++){
        cur_comb.push_back(i);
        combination(i,req_no,cur_comb);
        cur_comb.pop_back();
    }
}

vector<int> MyAlgo4::separation_oracle(int req_no, double &req_Us , vector<double>&alpha , map<pair<int,int>,double>&beta ){     
    vector<int> SPT;                  //nodes' parent in the spanning tree
    vector<int> best_path;
    double best_len=0; 
    int src = requests[req_no].get_source();
    int dst = requests[req_no].get_destination();


    SPT = Dijkstra(src, dst, req_no);                       //the first SPT is get by dijkstra
    int cur_node = src;                                     //counting the first path's U(X,Y)=c* e^r
    double c = 0;                                           //c=SUM[u,v]:alpha(u)+alpha(v)+beta(u,v)==X[u,v]
    double r = 0;                                           //r=SUM[u,v]:-ln(Pr(u,v))==Y[u,v]
    while(cur_node != dst){
        c += X(SPT[cur_node], cur_node, alpha, beta);  
        r += Y[req_no][{SPT[cur_node],cur_node}];  
        best_path.push_back(cur_node);
        cur_node = SPT[cur_node];
    } 
    best_path.push_back(dst);
    best_len = c * exp(r);                               

    map<pair<int, int>, bool> used_edge;
    vector<int> new_path;   
    pair<int,int> new_edge;

    for(unsigned int i = 0; i < SPT.size(); i++){
        int cur_node=i;
        while(cur_node!=dst){
            if(used_edge.find({cur_node,SPT[cur_node]})!=used_edge.end()){
                break;
            }
            used_edge[{cur_node,SPT[cur_node]}] = true;
            used_edge[{SPT[cur_node],cur_node}] = true;
            cur_node=SPT[cur_node];
        }
    }

    while(1){
        double minimum = numeric_limits<double>::infinity();
        for(int i = 0; i < graph.get_size(); i++){                 //creating many new SPT
            vector<int> neighbors = graph.get_neighbors_id(i);  
            for(auto neighbor : neighbors){
                double temp1 = 0, temp2 = 0;
                if(SPT[i] == neighbor || SPT[neighbor] == i){      // checking used edge or unused
                    continue;   
                }else{                                             // if unused
                    temp1 = X(i, neighbor, alpha , beta );
                    temp2 = Y[req_no][{i, neighbor}];
                    int cur_node = i;
                    while(cur_node != dst){
                        temp1 += X(cur_node, SPT[cur_node], alpha , beta );
                        temp2 += Y[req_no][{cur_node, SPT[cur_node]}];
                        cur_node = SPT[cur_node];
                    } 

                    cur_node = neighbor;
                    while(cur_node != dst){
                        temp1 -= X(cur_node, SPT[cur_node], alpha , beta );
                        temp2 -= Y[req_no][{cur_node, SPT[cur_node]}];
                        cur_node = SPT[cur_node];
                    }       
                    if(temp2 < 0 && temp1 > 0 ) {               // we need the smallest edge to change the SPT
                        if(i<neighbor){
                            if(used_edge.find({i, neighbor}) != used_edge.end()){
                                continue;
                            }
                        }
                        else{
                            if(used_edge.find({neighbor ,i}) != used_edge.end()){
                                continue;
                            }
                        }
                        if(minimum > - temp1 / temp2){
                            new_edge = {i, neighbor};
                            minimum = - temp1 / temp2;
                        }
                    }
                }
            }
        }        // 找到最小的 edge 

        if(minimum == numeric_limits<double>::infinity()){   //原本設計是有break,但之後用不到
            break;
        }else{
            new_path.clear();
        }
        used_edge[new_edge] = true;
        used_edge[{new_edge.second,new_edge.first}] = true;
        change_edge_num++;
        SPT[new_edge.second]=new_edge.first;

        cur_node = src;                                   
        while(cur_node != dst) {
            new_path.push_back(cur_node);
            cur_node = SPT[cur_node];
        }       
        new_path.push_back(dst);
        
        double new_len = 0;                                         //counting the new path's U(X,Y)=c* e^r
        c = 0;
        r = 0;
        for(unsigned int i = 0; i < new_path.size() - 1; i++){
            if(new_path[i] < new_path[i+1]){                        //[can alter]no need if else
                c += X(new_path[i], new_path[i+1], alpha , beta );
                r += Y[req_no][{new_path[i], new_path[i+1]}];
            }
            else{
                c += X(new_path[i+1], new_path[i], alpha , beta );  
                r += Y[req_no][{new_path[i+1], new_path[i]}];           
            }
            //cout<<"PATH:["<<new_path[i]<<" "<<new_path[i+1]<<"] with tot "<< c <<" / " << r <<endl;
        }

        new_len =  c * exp(r);
        if(new_len < best_len){
            best_len = new_len;
            best_path = new_path;                                            //路線修改,新的spt產生
        } 

    }
    req_Us += best_len;                                                    //we need k times this_time_U,add to req_Us
    return best_path;                                        
}

void MyAlgo4::find_bottleneck(vector<int> path, int req_no, bool obj_change, vector<double>&alpha , map<pair<int,int>,double>&beta){
    double min_s_u = numeric_limits<double>::infinity();
    double min_s_uv = numeric_limits<double>::infinity();
    vector<double> s_u(graph.get_size() + 5);
    vector<double> s_uv(graph.get_size() + 5);                                               

   
    for(unsigned int i = 0; i < path.size(); i++){
        if(i == 0 || i == path.size() - 1)
            s_u[path[i]] = graph.Node_id2ptr(path[i])->get_memory_cnt();        //if src or dst,then only cost 1
        else
            s_u[path[i]] = graph.Node_id2ptr(path[i])->get_memory_cnt() / 2.0;    //else cost 2
        if(s_u[path[i]] < min_s_u)
            min_s_u = s_u[path[i]];
    }

    for(unsigned int i = 0; i < path.size() - 1; i++){
        s_uv[i] = graph.get_channel_size(path[i], path[i+1]);                   //channel min
        if(s_uv[i] < min_s_uv)
            min_s_uv = s_uv[i];
    }

    int rate = 1;
    double s = min(min_s_u,min_s_uv);
    for(int i = 0; i < rate; i++){
        if(obj_change==true){
            if(x_i_p.find(path) != x_i_p.end())                                         //add flow to path
                x_i_p[path] += s;
            else
                x_i_p[path] = s;
        
            for(auto id : path){
                obj += (alpha[id] * (1 + epsilon * s / s_u[id]) - alpha[id]) * graph.Node_id2ptr(id)->get_memory_cnt();;            //alter alpha,beta
                alpha[id] = alpha[id] * (1 + epsilon * s / s_u[id]);
            }

            for(unsigned int i = 0; i < path.size() - 1; i++){
                obj += (beta[{path[i], path[i+1]}] * (1 + epsilon * s / s_uv[i]) -  beta[{path[i], path[i+1]}]) * graph.get_channel_size(path[i], path[i+1]);;
                beta[{path[i], path[i+1]}] = beta[{path[i], path[i+1]}] * (1 + epsilon * s / s_uv[i]);
                beta[{path[i+1], path[i]}] = beta[{path[i+1], path[i]}] * (1 + epsilon * s / s_uv[i]);
            }
       }
       else{
            for(auto id : path){
                alpha[id] = alpha[id] * (1 + epsilon * s / s_u[id]);
            }

            for(unsigned int i = 0; i < path.size() - 1; i++){
                beta[{path[i], path[i+1]}] = beta[{path[i], path[i+1]}] * (1 + epsilon * s / s_uv[i]);
                beta[{path[i+1], path[i]}] = beta[{path[i+1], path[i]}] * (1 + epsilon * s / s_uv[i]);
            } 
       }
    }
}

double MyAlgo4::changing_obj(){
    double temp_obj = 0.0;
    for(unsigned int i = 0; i < alpha.size(); i++){
        temp_obj += alpha[i] * graph.Node_id2ptr(i)->get_memory_cnt();
    }
    
    for(auto it : beta){
        temp_obj += it.second * graph.get_channel_size(it.first.first, it.first.second);
    }

    return temp_obj;
}

void MyAlgo4::find_violate(){
    vector<double> used_memory(graph.get_size());
    map<vector<int>, double> used_channel;
    map<pair<int, int>, double> used_request;

    for(auto &it : x_i_p){
        vector<int> path = it.first;
        int src = path[0];
        int dst = path.back();
    
        if(used_request.find({src, dst}) != used_request.end())
            used_request[{src, dst}] += it.second;
        else
            used_request[{src, dst}] = it.second;
        
        for(unsigned int i = 0; i < path.size() - 1; i++){
            used_memory[path[i]] += it.second;                         //memory add
            used_memory[path[i+1]] += it.second;
            if(path[i] < path[i+1]){
                auto iter = used_channel.find({path[i], path[i+1]});
                if(iter != used_channel.end()){    //channel add
                    used_channel[{path[i], path[i+1]}] += it.second;
                }
                else{
                    used_channel[{path[i], path[i+1]}] = it.second;
                }
            }
            else{
                auto iter = used_channel.find({path[i+1], path[i]});
                if(iter != used_channel.end()){
                    used_channel[{path[i+1], path[i]}] += it.second;
                }
                else{
                    used_channel[{path[i+1], path[i]}] = it.second;
                }  
            }
        }
    }


    double max_magni = 0.0;
    double cur_magni;

    for(auto it : used_request){
        int src = it.first.first;
        int dst = it.first.second;
        int req_no = -1;
        for(unsigned int i = 0; i < requests.size();i ++){
            if(requests[i].get_source() == src && requests[i].get_destination() == dst){
                req_no = i;
                break;
            }
        }
        cur_magni = it.second / requests[req_no].get_send_demand();
        if(cur_magni > max_magni){
            max_magni = cur_magni;
        }
    }

    for(auto it : used_channel){
        cur_magni = it.second / graph.get_channel_size(it.first[0],it.first[1]);
        if(cur_magni > max_magni){
            max_magni = cur_magni;
        }
    }

    for(int i = 0; i < graph.get_size(); i++){
        cur_magni = used_memory[i] / graph.Node_id2ptr(i)->get_memory_cnt();
        if(cur_magni > max_magni){
            max_magni = cur_magni;
        }
    }


    

    //cout << "Magnification:" << max_magni << endl;

    for(auto &x : x_i_p){

        x.second /= max_magni;
        
    }
    //check memory_and_channel
    /*
    for(unsigned int i=0;i<used_memory.size();i++){
        cout<<i<<" with memory "<<used_memory[i]<<endl;
    }
    for(auto it:used_channel){
        cout<<"["<<it.first[0]<<","<<it.first[1]<<"]:"<<it.second<<endl;
    }
    */
}


vector<map<vector<int>, int>> MyAlgo4::rounding(){
    vector<map<vector<int>, double>> each_request(requests.size());
    vector<map<vector<int>, int>> I_request(requests.size());
    for(auto it : x_i_p){
        vector<int> path = it.first;
        int src = path[0];
        int dst = path.back();
        for(unsigned int i = 0; i < requests.size(); i++){
            if(src == requests[i].get_source() && dst == requests[i].get_destination()){
                each_request[i][path] = it.second;
                break;
            }
        }
    }

    for(unsigned int i = 0; i < requests.size(); i++){
        double total_prob = 0;
        double used_prob = 0;
        int used_I=0;
        int distri_I=0;
        vector<double>accumulate;
        accumulate.push_back(0.0);                                              // [0,a1,a2,...,0]
        for(auto it : each_request[i]){                    
            double frac_prob;

            int i_prob = it.second;                                             //每個path先取整數部分=>確定分配
            I_request[i][it.first] = i_prob;
            used_I+=i_prob;

            frac_prob = it.second - i_prob;                                     //total_prob代表random區間,丟進accumulate
            total_prob += frac_prob;
            accumulate.push_back(total_prob);
            used_prob += it.second;
        }
        used_I += (int)(requests[i].get_send_demand()- used_prob);               //unused_I=取底[ri - sum(request.I) - (unused.I)]
        distri_I=requests[i].get_send_demand()-used_I;

        accumulate.push_back(0.0);
        // cout<<"total_prob:"<<total_prob<<" distri_I:"<<distri_I<<endl;
        // cout<<"accumulate:";
        // for(auto it:accumulate){
        //     cout<<it<<" ";
        // }
        // cout<<endl;
        for(int j = 0; j < distri_I; j++){
            random_device rd;  
            mt19937 gen(rd()); 
            uniform_real_distribution<double> dis(0.0, total_prob);
            double temp = dis(gen);
            for(unsigned int k = 0; k < accumulate.size() - 1; k++){
                // cout<<"ramdom:"<<temp<<endl;
                if(temp > accumulate[k] && temp < accumulate[k+1]){
                    unsigned int index = 0;
                    for(auto it : each_request[i]){
                        if(index == k){
                            I_request[i][it.first]++;
                            //cout<<"random prob="<<temp<<" so "<<k<<" added!\n";
                            break;
                        }
                        index += 1;
                    }
                    break;
                }
            }
        }
    }

    return I_request;
}

void MyAlgo4::check_enough(vector<map<vector<int>, int>> &path){
    vector<int> memory_used(graph.get_size());
    map<vector<int>,int> channel_used; 
    vector<int> over_memory(graph.get_size());
    map<vector<int>,int> over_channel;
    map<vector<int>,int>::iterator iter;
    for(int i = 0; i <(int)path.size(); i++){
        for(auto it : path[i]){
            vector<int> cur_path = it.first;
            for(int j = 0; j < (int)cur_path.size() - 1; j++){
                memory_used[cur_path[j]] += it.second;
                memory_used[cur_path[j+1]] += it.second;
                iter = channel_used.find({cur_path[j],cur_path[j+1]});
                if(iter != channel_used.end()){
                    channel_used[{cur_path[j], cur_path[j+1]}] += it.second;
                    channel_used[{cur_path[j+1], cur_path[j]}] += it.second;
                }
                else{
                    channel_used[{cur_path[j], cur_path[j+1]}] = it.second;
                    channel_used[{cur_path[j+1], cur_path[j]}] = it.second;
                }
            }
        }
    }

    for(int i = 0; i < graph.get_size(); i++){
        over_memory[i] = memory_used[i] - graph.Node_id2ptr(i)->get_memory_cnt();
        for(auto it : graph.get_neighbors_id(i)){
            iter = over_channel.find({i, it});
            if(iter != over_channel.end()){
               over_channel[{i, it}] -= graph.get_channel_size(i, it) / 2;
               over_channel[{it, i}] -= graph.get_channel_size(i, it) / 2; 
            }
            else{
               over_channel[{i, it}] = -graph.get_channel_size(i, it) / 2;
               over_channel[{it, i}] = -graph.get_channel_size(i, it) / 2;
            }
        }
    }

    for(auto &it : over_channel){
        iter = channel_used.find(it.first);
        if(iter != channel_used.end()){
            it.second = channel_used[{it.first}] + it.second ; 
        }
    }

    bool flag;
    while(1){
        flag = true;
        for(int i = 0; i < (int)over_memory.size(); i++){
            if(over_memory[i] > 0){
                flag = false;
            }
        }
        for(auto it : over_channel){
            if(it.second > 0){
                flag = false;
            }
        }
        if(flag == true){
            // cout << "before" << endl;
            // for(unsigned int i = 0; i < path.size(); i++){
            //     for(auto it:path[i]){
            //         vector<int>Final_path =it.first;
            //         for(auto it2:Final_path){
            //             cout<<it2<<" ";
            //         }
            //         cout<<"     Qubits:"<<it.second<<endl;
            //         requests[i].add_cur(it.second);
            //     }
            // }
            cout<<"--------------Reduece finish-------------\n";
            readd(path,over_memory,over_channel);  
            break;
        }
        int long_len = 0;
        int long_req = -1;
        vector<int> long_path;
        for(int i = 0; i < (int)path.size(); i++){
            for(auto it : path[i]){
                int associate_flag=false;
                /*
                for(auto temp:it.first){
                    cout<<temp<<" ";
                }
                cout<<"-----------"<<endl;
                */
                for(int j=0;j<(int)it.first.size()-1;j++){

                    //cout<<"memory check:"<<j<<"||"<<over_memory[it.first[j]]<<endl;
                    if(over_memory[it.first[j]]>0){
                        associate_flag=true;
                        break;
                    }
                    //cout<<"channel check:"<<j<<"/"<<j+1<<"||"<<over_channel[{it.first[j],it.first[j+1]}]<<endl;
                    iter = over_channel.find({it.first[j],it.first[j+1]});
                    if(iter!=over_channel.end() && over_channel[{it.first[j],it.first[j+1]}]>0){
                        associate_flag=true;
                        break;
                    }

                }
                if(over_memory[it.first[it.first.size()-1]]>0){
                    associate_flag=true;
                }

                if(associate_flag==true && (int)it.first.size() > long_len && it.second > 0){
                    long_len = it.first.size();
                    long_path = it.first;
                    long_req = i;
                }
            }
        }
        for(int i = 0; i < (int)long_path.size() - 1; i++){
            over_memory[long_path[i]]--;
            over_memory[long_path[i+1]]--;
            over_channel[{long_path[i], long_path[i+1]}]--;
            over_channel[{long_path[i+1], long_path[i]}]--;
        }
        path[long_req][long_path]--;
    }  
}  

void MyAlgo4::readd(vector<map<vector<int>, int>> &path,vector<int> &over_memory,map<vector<int>,int> &over_channel){
    for(unsigned int i = 0; i < path.size(); i++){
        for(auto it : path[i]){
            requests[i].add_cur(it.second);
        }
    }
    vector<pair<vector<int>, int>> re;
    int max = -1;
    for(unsigned int i = 0; i < requests.size(); i++){
        for(auto it : path[i]){
            if(max < it.second){
                max = it.second;
            }
        }
    }

    for(int i = max; i >= 0; i--){
        for(unsigned int j = 0; j < requests.size(); j++){
            for(auto it : path[j]){
                if(i == it.second){
                    re.push_back({it.first, j});
                }
            }
        }
    }
    bool flag = true;
    while(flag){
        flag = false;
        for(unsigned int i = 0; i < re.size(); i++){
            if(requests[re[i].second].get_send_demand() > requests[re[i].second].get_cur_send()){
                vector<int> each_path = re[i].first;
                bool assign = true;
                for(unsigned int j = 0; j < each_path.size() - 1; j++){
                    if(j == 0){
                        if(over_memory[each_path[j]] >= 0){
                            assign = false;
                        }
                    }
                    else{
                        if(over_memory[each_path[j]] >= -1){
                            assign = false;
                        }
                    }
                    if(over_channel[{each_path[j],each_path[j+1]}] >= 0){
                        assign = false;
                    }
                }
                if(over_memory[each_path[each_path.size()-1]] >= 0){
                    assign = false;
                }
                if(assign == true ){
                    requests[re[i].second].add_cur(1);
                    for(auto it : path[re[i].second]){
                        if(it.first == re[i].first){
                            path[re[i].second][it.first] += 1;
                            cout << "!!PATH +++" << endl;
                            flag = true;
                            for(int j = 0; j < (int)each_path.size() - 1; j++){
                                over_memory[each_path[j]]++;
                                over_memory[each_path[j+1]]++;
                                over_channel[{each_path[j], each_path[j+1]}]++;
                                over_channel[{each_path[j+1], each_path[j]}]++;
                            }
                            break;
                        }
                    }
                }
            }
        }
    }
    // for(auto x : over_memory){
    //     cout << "node: " << x << endl;
    // }
    // for(auto x : over_channel){
    //     cout<< "EDGE: ";
    //     for(auto a : x.first ){
    //         cout<< a << " ";
    //     }
    //     cout << x.second << endl;
    // }

}

void MyAlgo4::dfs(int src, int dst, vector<vector<int>> &ans, vector<int> &path, vector<bool> &visited ){
        //base case
    visited[src] = true;
    path.push_back(src);
    if(src == dst){
        ans.push_back(path);
    } 

    else{
        for(auto i : graph.get_neighbors_id(src)){ 
            if(!visited[i]){
                dfs(i, dst, ans, path, visited);
            }
        }
    }
    visited[src] = false;
    path.pop_back();

}


void MyAlgo4::calculate(){
    double sum=0.0;

    int t = 1;

    for(auto it:x_i_p){
        double prob=1;
        vector<int>path=it.first;
        for(unsigned int i=0;i < it.first.size() - 1;i++){
            prob*=exp(graph.Node_id2ptr(path[i])->distance(*graph.Node_id2ptr(path[i+1]))*(-graph.get_entangle_alpha()));
        }
        for(unsigned int i=1;i<it.first.size()-1;i++){
            prob*=graph.Node_id2ptr(path[i])->get_swap_prob();  
        }
        sum+=it.second*prob;

        t++;
    }
    cerr << "sum = " << sum << endl;
    res["primal"] = sum / (1 - epsilon) / (1 - epsilon);
}

vector<vector<int>> MyAlgo4::allPathsSourceTarget(int src, int dst){
    vector<vector<int>> ans;
    vector<int> path;
    vector<bool> visited(graph.get_size());
    for(int i = 0; i < graph.get_size(); i++){
        visited[i] = false;
    }
    dfs(src, dst, ans, path, visited);
    return ans;
}

vector<map<vector<int>, int>> MyAlgo4::Greedy_rounding(){
    vector<map<vector<int>, double>> each_request(requests.size());
    vector<map<vector<int>, int>> I_request(requests.size());
    for(auto it : x_i_p){
        vector<int> path = it.first;
        int src = path[0];
        int dst = path.back();
        for(unsigned int i = 0; i < requests.size(); i++){
            if(src == requests[i].get_source() && dst == requests[i].get_destination()){
                each_request[i][path] = it.second;
                break;
            }
        }
    }

    for(unsigned int i = 0; i < each_request.size(); i++){
        for(auto it:each_request[i]){
            vector<int>undistri_path =it.first;
        }
    }
	vector<int> used_I(requests.size());										//第 i 個 request 目前用了幾調 path
	vector< tuple<double, int, vector<int>> > fractional_xip;	
    for(unsigned int i = 0; i < requests.size(); i++){
        used_I[i] = 0;
		for(auto it : each_request[i]){                    
            double frac_prob;

            int i_prob = it.second;                                             //每個path先取整數部分=>確定分配
            I_request[i][it.first] = i_prob;
            used_I[i] += i_prob;
			assign_resource(it.first, i_prob, i);
            frac_prob = it.second - i_prob;                                     //total_prob代表random區間,丟進accumulate
			fractional_xip.emplace_back(frac_prob*graph.find_success_probability(it.first), i, it.first);
        }                                                   //unused_I=取底[ri - sum(request.I) - (unused.I)]
    }

	// 對 x^i_p 由小到大排序
	sort(fractional_xip.begin(), fractional_xip.end());
	reverse(fractional_xip.begin(), fractional_xip.end());

	// 如果資源足夠讓 x 變成 1 ，就直接讓 x 變成 1 
	for(auto it:fractional_xip){
		vector<int> extra_path;
		double x_prob;
		int request_id;
		tie(x_prob, request_id, extra_path) = it;
		if(find_width(extra_path) >= 1 && used_I[request_id] < requests[request_id].get_send_demand()){
			assign_resource(extra_path, 1, request_id);
			used_I[request_id] += 1;
            I_request[request_id][extra_path]++;
		}
	}
	// 如果還有剩下資源的話，盡量塞爆
	for(auto it:fractional_xip){
		vector<int> extra_path;
		double x_prob;
		int request_id;
		tie(x_prob, request_id, extra_path) = it;
		int width = 0;
		int extra_send_demand = requests[request_id].get_send_demand() - used_I[request_id];
		width = min(find_width(extra_path), extra_send_demand);
		if(width >= 1){
			assign_resource(extra_path, width, request_id);
            used_I[request_id] += width;
			I_request[request_id][extra_path]++;
		}
	}
	for(int request_id=0;request_id<(int)requests.size();request_id++){
		// cerr<<"GG: It work?"<<endl;
		while(requests[request_id].get_send_demand() - used_I[request_id] > 0){
			vector<int> extra_path = BFS(requests[request_id].get_source(), requests[request_id].get_destination());
			int width = 0;
			if(extra_path.size() != 0){
				width = min(find_width(extra_path), requests[request_id].get_send_demand() - used_I[request_id]);
				assign_resource(extra_path, width, request_id);
				used_I[request_id] += width;
			}
			if(width == 0){
				break;
			}
		}
	}
    return I_request;

}

void MyAlgo4::path_assignment(){
    for(int i = 0; i < graph.get_size(); i++){
        int mem = graph.Node_id2ptr(i)->get_memory_cnt();
        graph.Node_id2ptr(i)->revise(mem);
    }
    
    initialize();
    
    MyAlgo4::comb.resize(requests.size());
    
    for(unsigned int i = 0; i < requests.size(); i++){
        int src = requests[i].get_source();
        int dst = requests[i].get_destination();
        all_source_target_path.push_back(allPathsSourceTarget(src, dst));
        vector<int> cur_comb;
        combination(0,i,cur_comb);
    }
    
    obj = M * delta;
    vector<vector<int>> best_path;

    while(obj < 1){
        int req_no = 0;
        double smallest_U = numeric_limits<double>::infinity();       //this compare all request U
        vector<double> U;
        vector<vector<vector<int>>>all_path;
        //vector<vector<vector<int>>>all_path2;                      //compare if finding k path = opt
        all_path.resize(requests.size());
        //all_path2.resize(requests.size());
        U.resize(requests.size());
        // cout<<"\n------New round-------\n";
        
        /*     [FOR 暴力找最佳SET]
        //#pragma omp parallel for
        for(unsigned int i = 0; i < requests.size(); i++){
            vector<int>best_comp;
            double min_U= numeric_limits<double>::infinity();         //this compare one request U <temporary use>
            for(unsigned int j = 0; j < comb[i].size(); j++){
                vector<double>alpha_temp=alpha;
                map<pair<int,int>,double>beta_temp=beta;
                double c=0;
                double r=0;
                //cout<<"CombInfo"<<j<<endl;
                for(unsigned int k = 0;k<comb[i][j].size(); k++){
                    vector<int> cur_path=MyAlgo4::all_source_target_path[req_no][comb[i][j][k]];
                    //cout<<"PathInfo"<<k<<":";
                    for(unsigned int l=0;l<cur_path.size()-1;l++){
                        //cout<<cur_path[l]<<" ";
                        c += X(cur_path[l], cur_path[l+1], alpha_temp , beta_temp );  
                        r += Y[req_no][{cur_path[l], cur_path[l+1]}]; 
                    }
                    //cout<<cur_path.back()<<endl;
                    find_bottleneck(cur_path, i, false, alpha_temp, beta_temp);
                    //cout<<"process2"<<endl;
                    //cout<<c<<" "<<r<<endl;
                }
                
                //cout<<"CombValue:"<<c*exp(r)<<endl<<endl;
                if(c * exp(r) < min_U){
                    min_U = c * exp(r);
                    best_comp=comb[i][j]; 
                    //cout<<"Request"<<i<<"-Comb"<<j<<":"<<c*exp(r)<<endl;
                }
                //cout<<"process3"<<endl;
            }
            for(auto  it:best_comp){
                all_path2[i].push_back(MyAlgo4::all_source_target_path[i][it]);
            }
        }
        */

        #pragma omp parallel for                                          //這樣會不會互相影響到不同的request?
        for(unsigned int i = 0; i < requests.size(); i++){                //for each request,do their k separation oracle          
            vector<double>alpha_temp=alpha;                              
            map<pair<int,int>,double>beta_temp=beta;

            for(int j = 0; j < 3/*this is k*/; j++){
                all_path[i].push_back(separation_oracle(i, U[i],alpha_temp,beta_temp));
                find_bottleneck(all_path[i].back(), i, false, alpha_temp, beta_temp);
            }
        }

        for(unsigned int i = 0; i < requests.size(); i++){
            if(U[i] < smallest_U){
                smallest_U  = U[i];
                best_path = all_path[i];
                req_no = i;
            }
        } 
        /*
        for(auto it:all_path){
            for(auto it2:it){
                cout<<"[1]";
                for(auto it3:it2){
                    cout<<it3<<" ";
                }
                cout<<endl;
            }
        }
        */

        /*     [FOR 暴力找最佳SET]
        for(auto it:all_path2){
            for(auto it2:it){
                cout<<"[2]";
                for(auto it3:it2){
                    cout<<it3<<" ";
                }
                cout<<endl;
            }
        }
        */

        for(auto it:all_path[req_no]){
            find_bottleneck(it, req_no, true, MyAlgo4::alpha, MyAlgo4::beta);
        }
        
        cout << obj << endl;
    }

    // calculate();
    find_violate();
    calculate();
    vector<map<vector<int>, int>>path = Greedy_rounding();
    res["change_edge_num"] = change_edge_num;
    res["diff_edge_num"] = diff_num;

}   

