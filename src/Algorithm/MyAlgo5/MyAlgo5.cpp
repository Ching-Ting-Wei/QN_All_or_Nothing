#include "MyAlgo5.h"

MyAlgo5::MyAlgo5(string filename, int request_time_limit, int node_time_limit, double swap_prob, double entangle_alpha)
    :AlgorithmBase(filename, "MyAlgo5", request_time_limit, node_time_limit, swap_prob, entangle_alpha , true /*let it be alway true*/){
    if(DEBUG) cerr<<"new MyAlgo5"<<endl;
}
MyAlgo5::MyAlgo5(string filename, int request_time_limit, int node_time_limit, double swap_prob, double entangle_alpha, double epsilon)
    :AlgorithmBase(filename, "MyAlgo5_" + to_string(epsilon) , request_time_limit, node_time_limit, swap_prob, entangle_alpha , true /*let it be alway true*/), epsilon(epsilon){
    if(DEBUG) cerr<<"new MyAlgo5"<<endl;
}


MyAlgo5::~MyAlgo5(){
    if(DEBUG) cerr << "delete MyAlgo5" << endl;
}

void MyAlgo5::entangle(){
    AlgorithmBase::base_entangle();
}

void MyAlgo5::swap(){
     AlgorithmBase::base_swap();
}

void MyAlgo5::send(){
     AlgorithmBase::base_send();
}

void MyAlgo5::next_time_slot(){
     AlgorithmBase::base_next_time_slot();
}

double MyAlgo5::X(int u, int v, int i){
	double weight = alpha[u] + alpha[v] + beta[{u, v}];
	if(requests[i].get_source() == u || requests[i].get_source() == v){
		weight += tau[i];
	}
	return weight;
}

void MyAlgo5::yen(int src,int dst,int K,int req_no){
    vector<vector<int>> copy_graph;
    vector<vector<int>> candidate_path;
    vector<int> parent;
    vector<int> rootpath;
    vector<int> spurpath;
    vector<int> totalpath;
    vector<pair<int,int>> get_delete;
    int spurnode;
    double min_dist=numeric_limits<double>::infinity();
    double cur_dist;
    for(int i = 0 ;i < graph.get_size() ;i++ ){                                                     //create graph to alter
        copy_graph.push_back(graph.get_neighbors_id(i));  
    }

    all_given_path[req_no].push_back(Dijkstra(copy_graph,src,dst,req_no,rootpath));                 //get the first path
    /*
    cout<<"[First path]";
    for(auto it:all_given_path[req_no][0]){
        cout<<it<<" ";
    }
    cout<<endl;
    */
    for(int k = 1 ;k < K ;k++ ){                                                                    //find k path
        //cout<<"Noe finding the "<< k <<"path............."<<endl;
        for(int i = 0; i <= (int)(all_given_path[req_no][k - 1].size()) - 2 ; i++){                 //spurnode index     
            //cout<<"[i]:"<<i<<" [all_given.size]"<<"["<<k<<"]"<<all_given_path[req_no][k-1].size()<<endl;
            spurnode = all_given_path[req_no][k - 1][i];
            for(int j = 0; j <= i; j++){
                rootpath.push_back(all_given_path[req_no][k-1][j]);                                 //rootpath done
            }
                
            for(auto it:all_given_path[req_no]){                                                    //if 重疊，在graph刪除(i,i+1)
                int index=0;
                bool same=true;
                //cout<<"rootpath:";
                for(auto it2:rootpath){
                    //cout<<it2<<" ";
                    if(it2!=it[index++]){
                        same=false;
                    }
                }
                //cout<<endl;
                if(same){                                                                           //overlap,then remove edge (i,i+1)
                    vector<int>::iterator iter = find(copy_graph[it[i]].begin(), copy_graph[it[i]].end(), it[i+1]);        
                    //cout<<"[What iter find]:"<<*iter<<endl;
                    if(iter != copy_graph[it[i]].end()){
                       copy_graph[it[i]].erase(iter); 
                    }
                    
                    iter = find(copy_graph[it[i+1]].begin(), copy_graph[it[i+1]].end(), it[i]);
                    if(iter != copy_graph[it[i+1]].end()){
                        copy_graph[it[i+1]].erase(iter);
                        get_delete.push_back({it[i],it[i+1]});
                    }
                    
                }      
            }
            spurpath = Dijkstra(copy_graph,spurnode, dst, req_no, rootpath);                         //rootpathnodes are banned 
            totalpath=rootpath;                                                                     //totalpath = rootpath + spurpath
            totalpath.pop_back();
            rootpath.clear();
            //cout<<"spurpath";
            for(auto it:spurpath){
                //cout<<it<<" ";
                totalpath.push_back(it);
            }
            //cout<<endl;

            bool inB=false;                                                                         //check if totalpath is already in candidate path
            for(auto it:candidate_path){
                if(it == totalpath){
                    inB=true;
                }
            }
            if(!inB){
                candidate_path.push_back(totalpath);
            }
                
            for(auto it:get_delete){                                                                //restore edges
                copy_graph[it.first].push_back(it.second);
                copy_graph[it.second].push_back(it.first);
            }
            get_delete.clear();
        }        
        // for(auto it:candidate_path){
        //     cout<<"Candidate_path";
        //     for(auto it2:it){
        //         cout<<it2<<" ";
        //     }
        //     cout<<endl;
        // }

        if(candidate_path.size() == 0){
            break;
        }
        
        for(auto it:candidate_path){                                                                //find the shortest path in candidate path
            cur_dist=0;
            for(unsigned int i = 0; i < it.size()-1; i++){
                cur_dist += graph.Node_id2ptr(it[i])->distance(*graph.Node_id2ptr(it[i+1]));
            }
            if(cur_dist<min_dist){
                min_dist=cur_dist;
                totalpath=it;
            }
        }
        all_given_path[req_no].push_back(totalpath);
        /*
        cout<<"[PATH]";
        for(auto it:totalpath){
            cout<<it<<" ";
        }
        cout<<endl;
        */
        for(auto it:candidate_path){
            it.clear();
        }
        
        candidate_path.clear();
    }    
}

void MyAlgo5::initialize(){                                               //Main idea [Path is edge]

    M = graph.get_size() + graph.get_num_of_edge() + requests.size();              //M=V+E+I
    // delta = (1 + epsilon) * pow(((1 + epsilon) * M), (-1 / epsilon));         
    delta = (1 + epsilon) * (1.0 / pow((1 + epsilon) * M, 1.0 / epsilon));

    for(int i = 0; i < graph.get_size(); i++){
        alpha.emplace_back(delta / graph.Node_id2ptr(i)->get_memory_cnt());       //alpha=delta/Mu

        vector<int> temp = graph.get_neighbors_id(i);                             //beta=delta/Cuv
        for(auto it: temp){
            beta[make_pair(i,it)] = delta / (graph.get_channel_size(i, it));
            beta[make_pair(it,i)] = delta / (graph.get_channel_size(i, it));
        }
    }

    for(unsigned  i = 0; i < requests.size(); i++){                               //tau=delta/ri
        tau.emplace_back(delta);
    }
    
    Y.resize(requests.size() + 5);
    for(int i = 0; i < graph.get_size(); i++){
        vector<int> temp = graph.get_neighbors_id(i);                             
        for(auto it: temp){

            for(unsigned  j = 0; j < requests.size(); j++){                       //Y=-ln(edge)-ln(repeater_1^(1/2))-ln(repeater_2^(1/2))
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

vector<int> MyAlgo5::Dijkstra(vector<vector<int>>&copy_graph, int src, int dst, int req_no, vector<int>&get_delete){ 
    const double INF = numeric_limits<double>::infinity();
    int n = copy_graph.size();
    vector<double> distance(n, INF);
    vector<int> parent(n, -1);
    vector<bool> used(n, false);
    priority_queue<pair<double, int>, vector<pair<double, int>>, greater<pair<double, int>>> pq;
    
    distance[src] = 0;
    pq.push({0, src});
    
    while(!pq.empty()) {
        int cur_node = pq.top().second;
        pq.pop();
        if(used[cur_node]) continue;
        used[cur_node] = true;
        for(auto neighbor : copy_graph[cur_node]){
            bool jump_out = false;
            for(unsigned int i = 0; i < get_delete.size(); i++){
                if(neighbor == get_delete[i] && i != get_delete.back()){
                    jump_out = true;
                }
            }
            if(jump_out){
                continue;
            }
            if(distance[cur_node] + graph.Node_id2ptr(cur_node)->distance(*graph.Node_id2ptr(neighbor)) < distance[neighbor]) {
                distance[neighbor] = distance[cur_node] + graph.Node_id2ptr(cur_node)->distance(*graph.Node_id2ptr(neighbor));
                parent[neighbor] = cur_node;
                pq.push({distance[neighbor], neighbor});
            }
        }
    }
    
    if(distance[dst] >= INF){
        cout<<"[Dijkstra Not Find Route]"<<endl;
        return{};
    }
    int cur_node = dst;
    vector<int> path;
    while(cur_node != -1) {
        path.push_back(cur_node);
        cur_node = parent[cur_node];
    }
    reverse(path.begin(), path.end());
    return path;
    
}       
   
vector<int> MyAlgo5::SepDijkstra(int src, int dst, int req_no, vector<vector<double>> &path_graph_X){ 

    int path_num = 10;
    vector<vector<double>> transpose_graph(qubit_num * path_num + 2, vector<double>(qubit_num * path_num  + 2));
    for (int i = 0; i < qubit_num * path_num + 2; i++)
        for (int j = 0; j < qubit_num * path_num + 2; j++)
            transpose_graph[i][j] = path_graph_X[j][i];
    


    const double INF = numeric_limits<double>::infinity();
    vector<double> distance(qubit_num * path_num + 2, INF);
    vector<int> parent(qubit_num * path_num + 2, -1);
    vector<bool> used(qubit_num * path_num + 2, false);
    priority_queue<pair<double, int>, vector<pair<double, int>>, greater<pair<double, int>>> pq;

    distance[dst] = 0;
    pq.push({0, dst});
    while(!pq.empty()) {
        int cur_node = pq.top().second;
        pq.pop();
        if(used[cur_node]) continue;
        used[cur_node] = true;
        for(int i=0; i < qubit_num * path_num + 2;i++) {
            if(transpose_graph[cur_node][i]!= -1){
                if(cur_node==21 && i == 31){
                    cout<<"[WHAT happend]"<< distance[cur_node] + transpose_graph[cur_node][i] <<endl;
                }
                if(distance[cur_node] + transpose_graph[cur_node][i] < distance[i]) {
                    distance[i] = distance[cur_node] + transpose_graph[cur_node][i];
                    parent[i] = cur_node;
                    
                    pq.push({distance[i], i});
                }
            }

        }
    }
    if(distance[dst] >= INF) return{};
    return parent;
}     

vector<int> MyAlgo5::separation_oracle(int req_no, double &req_Us, vector<vector<vector<double>>> &path_graph_X, vector<vector<vector<double>>> &path_graph_Y){  
   
    vector<int> SPT;                  //nodes' parent in the spanning tree
    vector<int> best_set;
    int path_num = 10;
    double best_len; 
    int src = 0;
    int dst = path_num * qubit_num + 1;

    SPT = SepDijkstra(src, dst, req_no, path_graph_X[req_no]);         //the first SPT is get by dijkstra
    cout<<"[SPT]:";
    for(auto it:SPT){
        cout<<it<<" ";
    }
    cout<<endl;
    int cur_node = src;                                     //counting the first path's U(X,Y)=c* e^r
    double c = 0;                                           //c = SUM[u,v]:alpha(u)+alpha(v)+beta(u,v)==X[u,v]
    double r = 0;                                           //r = SUM[u,v]:-ln(Pr(u,v))==Y[u,v]
    while(cur_node != dst){
        // if(cur_node < SPT[cur_node]){                       //[can alter]no need if else
        //     c += X[{cur_node,SPT[cur_node]}][req_no];               
        //     r += Y[req_no][{cur_node,SPT[cur_node]}];
        // }
        // else{
        //     c += X[{SPT[cur_node],cur_node}][req_no];  
        //     r += Y[req_no][{SPT[cur_node],cur_node}];           
        // }
        cout<<"[while]:"<<cur_node<<endl;
        c += path_graph_X[req_no][cur_node][SPT[cur_node]]; 
        r += path_graph_Y[req_no][cur_node][SPT[cur_node]];   
        best_set.push_back(cur_node);
        cur_node = SPT[cur_node];
    } 
    best_set.push_back(dst);
    best_len = c * exp(r);
    req_Us = best_len;

    // cout << "origin path: ";
    // for(auto p : best_set){
    //         cout << p << "->";
    // }
    // cout << endl;

    map<pair<int, int>, bool> used_edge;
    vector<int> new_path;   
    pair<int,int> new_edge;

    for(unsigned int i = 0; i < SPT.size(); i++){
        int cur_node = i;
        while(cur_node != dst){
            if(used_edge.find({cur_node,SPT[cur_node]}) != used_edge.end() || SPT[cur_node] == -1){
                break;
            }
            used_edge[{cur_node,SPT[cur_node]}] = true;
            // used_edge[{SPT[cur_node],cur_node}] = true;
            cur_node = SPT[cur_node];
        }
    }
    for(auto it:used_edge){
        cout<<it.first.first<<","<<it.first.second<<" is "<< it.second<<endl;
    }
    while(1){
        double minimum = numeric_limits<double>::infinity();
        for(int i = 0; i < path_num * qubit_num + 2; i++){                 //creating many new SPT
            for(int j = i + 1; j < path_num * qubit_num + 2; j++){
                double temp1 = 0, temp2 = 0;
                if(path_graph_X[req_no][i][j] != -1){
                    if(SPT[i] == j){      // checking used edge or unused
                        continue;   
                    }else{
                        temp1 = path_graph_X[req_no][i][j];
                        temp2 = path_graph_Y[req_no][i][j];
                        int cur_node = i;
                        while(cur_node != dst){
                            temp1 += path_graph_X[req_no][cur_node][SPT[cur_node]];
                            temp2 += path_graph_Y[req_no][cur_node][SPT[cur_node]];
                            cur_node = SPT[cur_node];
                        } 

                        cur_node = j;
                        while(cur_node != dst){
                            temp1 -= path_graph_X[req_no][cur_node][SPT[cur_node]];
                            temp2 -= path_graph_Y[req_no][cur_node][SPT[cur_node]];
                            cur_node = SPT[cur_node];
                        }
                        if(temp2 < 0 && temp1 > 0) {
                            if(used_edge.find({i,j}) != used_edge.end()){
                                continue;
                            }
                            if(minimum > -temp1 / temp2){
                                new_edge = {i,j};
                                minimum = -temp1 / temp2;
                            }
                        }

                    }   
                }
            }
        }        // 找到最小的 edge 
        
        cout<<"[new_edge]"<<new_edge.first<<","<<new_edge.second<<endl;
        if(minimum == numeric_limits<double>::infinity()){   //原本設計是有break,但之後用不到
            break;
        }else{
            new_path.clear();
        }
        used_edge[new_edge] = true;
        change_edge_num++;
        SPT[new_edge.first]=new_edge.second;    //***first to second?
 
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
            c += path_graph_X[req_no][cur_node][SPT[cur_node]]; 
            r += path_graph_Y[req_no][cur_node][SPT[cur_node]]; 
        }
        new_len =  c * exp(r);
        if(new_len < best_len){
            best_len = new_len;
            req_Us = best_len;
            best_set = new_path;                                            //路線修改,新的spt產生
        } 
        
    }    
    return best_set;                                                  
}

void MyAlgo5::find_bottleneck(vector<int> path, int req_no){
    
    double min_s_u = numeric_limits<double>::infinity();
    double min_s_uv = numeric_limits<double>::infinity();
    double s_i = requests[req_no].get_send_limit();                             //request[no] min
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

    int rate = 10;
    double s = min(min_s_u, min(min_s_uv, s_i));
    for(int i = 0; i < rate; i++){
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
        obj += (tau[req_no] * (1 +epsilon * s) - tau[req_no]);
        tau[req_no] = tau[req_no] * (1 + epsilon * s);
    }
}

double MyAlgo5::changing_obj(){
    double temp_obj = 0.0;
    for(unsigned int i = 0; i < alpha.size(); i++){
        temp_obj += alpha[i] * graph.Node_id2ptr(i)->get_memory_cnt();
    }
    
    for(auto it : beta){
        temp_obj += it.second * graph.get_channel_size(it.first.first, it.first.second);
    }

    for(unsigned int i = 0;i < requests.size(); i++){
        temp_obj += tau[i];
    }
    return temp_obj;
}

void MyAlgo5::find_violate(){
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
        cur_magni = it.second / requests[req_no].get_send_limit();
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

void MyAlgo5::check_enough(vector<map<vector<int>, int>> &path){
    vector<int> memory_used(graph.get_size());
    map<vector<int>,int> channel_used; 
    vector<int> over_memory(graph.get_size());
    map<vector<int>,int> over_channel;
    map<vector<int>,int>::iterator iter;
    for(int i = 0; i <(int)path.size(); i++){
        for(auto it : path[i]){
            vector<int> cur_set = it.first;
            for(int j = 0; j < (int)cur_set.size() - 1; j++){
                memory_used[cur_set[j]] += it.second;
                memory_used[cur_set[j+1]] += it.second;
                iter = channel_used.find({cur_set[j],cur_set[j+1]});
                if(iter != channel_used.end()){
                    channel_used[{cur_set[j], cur_set[j+1]}] += it.second;
                    channel_used[{cur_set[j+1], cur_set[j]}] += it.second;
                }
                else{
                    channel_used[{cur_set[j], cur_set[j+1]}] = it.second;
                    channel_used[{cur_set[j+1], cur_set[j]}] = it.second;
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

void MyAlgo5::readd(vector<map<vector<int>, int>> &path,vector<int> &over_memory,map<vector<int>,int> &over_channel){
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
            if(requests[re[i].second].get_send_limit() > requests[re[i].second].get_cur_send()){
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

void MyAlgo5::dfs(int src, int dst, vector<vector<int>> &ans, vector<int> &path, vector<bool> &visited){
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


void MyAlgo5::calculate(){
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

vector<map<vector<int>, int>> MyAlgo5::Greedy_rounding(){
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
		if(find_width(extra_path) >= 1 && used_I[request_id] < requests[request_id].get_send_limit()){
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
		int extra_send_limit = requests[request_id].get_send_limit() - used_I[request_id];
		width = min(find_width(extra_path), extra_send_limit);
		if(width >= 1){
			assign_resource(extra_path, width, request_id);
            used_I[request_id] += width;
			I_request[request_id][extra_path]++;
		}
	}
	for(int request_id=0;request_id<(int)requests.size();request_id++){
		// cerr<<"GG: It work?"<<endl;
		while(requests[request_id].get_send_limit() - used_I[request_id] > 0){
			vector<int> extra_path = BFS(requests[request_id].get_source(), requests[request_id].get_destination());
			int width = 0;
			if(extra_path.size() != 0){
				width = min(find_width(extra_path), requests[request_id].get_send_limit() - used_I[request_id]);
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

void MyAlgo5::create_pathGraph(vector<vector<vector<double>>> &path_graph_X, vector<vector<vector<double>>> &path_graph_Y, vector<vector<double>> &X_value, vector<vector<double>> & Y_value)
{
    // We need to build direct graph
    // [0] represent source node
    // [n*k+2] represent destination node
    int path_num = all_given_path[0].size();
    
    for(unsigned int i = 0; i < requests.size(); i++ ){
        for(unsigned int j = 0; j < all_given_path[i].size(); j++){
            double c = 0;
            double r = 0;
            for(unsigned int k = 0; k < all_given_path[i][j].size() - 1; k++){
                c += X(all_given_path[i][j][k], all_given_path[i][j][k+1], i);               
                r += Y[i][{all_given_path[i][j][k], all_given_path[i][j][k+1]}]; 
            }
            X_value[i][j] = c;
            Y_value[i][j] = r;
        }

    

        for(int j = 0; j < path_num - qubit_num + 1; j++){
            path_graph_X[i][0][j+1] = X_value[i][j];
            path_graph_Y[i][0][j+1] = Y_value[i][j];
        }

        for(int j = 0; j < qubit_num - 1; j++){
            for(int k = 0; k < path_num; k++){
                if((path_num - k) - (qubit_num - j) < 0) continue;
                for(int l = k + 1; l < path_num; l++){
                    if((path_num - l) - (qubit_num - 1 - j) < 0) continue;
                    path_graph_X[i][j * path_num + k + 1][(j + 1) * path_num + l + 1] = X_value[i][l];
                    path_graph_Y[i][j * path_num + k + 1][(j + 1) * path_num + l + 1] = Y_value[i][l];
                }
            }
        }   


        for(int j = 0; j < path_num; j++){
            path_graph_X[i][path_num * (qubit_num - 1) + j + 1][path_num * qubit_num + 1] = 0;
            path_graph_Y[i][path_num * (qubit_num - 1) + j + 1][path_num * qubit_num + 1] = 1;
        }
       
    }

 

    for(unsigned int i = 0; i < requests.size(); i++ ){
        for(int j = 0; j < path_num * qubit_num + 2; j++){
            for(int k = 0; k < path_num * qubit_num + 2; k++){
                if(path_graph_X[i][j][k] != -1){
                    cout << "1";
                }else{
                    cout << 0;
                }
            }
            cout << endl;
        }
        cout << "---------" << endl;
    }

}

void MyAlgo5::path_assignment(){

    all_given_path.resize(requests.size());
    for(unsigned int i = 0; i < requests.size(); i++){
        cout<<"[request]:"<<i<<endl;
        yen(requests[i].get_source(),requests[i].get_destination(),10,i);
    }

    /*
    show given path
    for(auto it:all_given_path){
        cout<<"Next Request -----------"<<endl;
        for(auto it2:it){
            cout<<"[Path]";
            for(auto it3:it2){
                cout<<it3<<" ";
            }
            cout<<endl;
        }
        cout<<endl;
    }
    */ 

   int path_num = 10;

    initialize();
    

    vector<vector<vector<double>>> path_graph_X(requests.size(), vector<vector<double>>(10 * qubit_num + 2, vector<double>(10 * qubit_num + 2,-1)));
    vector<vector<vector<double>>> path_graph_Y(requests.size(), vector<vector<double>>(10 * qubit_num + 2, vector<double>(10 * qubit_num + 2,-1)));
    vector<vector<double>> X_value(requests.size(), vector<double>(path_num));
    vector<vector<double>> Y_value(requests.size(), vector<double>(path_num));
    create_pathGraph(path_graph_X, path_graph_Y, X_value, Y_value);

    obj = M * delta;
    vector<int> best_set;
    vector<int> cur_set;

    while(obj < 1){
        int req_no = 0;
        double smallest_U = numeric_limits<double>::infinity();
        vector<double> U;
        vector<vector<int>>all_path_set;
        all_path_set.resize(requests.size());
        U.resize(requests.size());
        // cout<<"\n------New round-------\n";
        //#pragma omp parallel for
        for(unsigned int i = 0; i < requests.size(); i++){
            all_path_set[i] =  separation_oracle(i, U[i], path_graph_X, path_graph_Y);
            //cout << "smallest_U: " << smallest_U << " U: " << U << "\n\n"; 
        }
        for(auto it:all_path_set){
            cout<<"[SET]";
            for(auto it2:it){
                cout<<it2<<" ";
            }
            cout<<endl;
        }
        for(unsigned int i = 0; i < requests.size(); i++){
            //cout << "smallest_U: " << smallest_U << " U: " << U << "\n\n"; 
            if(U[i] < smallest_U){
                smallest_U  = U[i];
                best_set = all_path_set[i];
                req_no = i;
            }
        } 
        // cout << smallest_U << endl;

        find_bottleneck(best_set, req_no);
        
        cout << obj << endl;
        // obj = changing_obj();
        // cout<<"changing_obj obj: " << obj << endl ;
    }
/*
    // calculate();
    find_violate();
    calculate();
    vector<map<vector<int>, int>>path = Greedy_rounding();
    res["change_edge_num"] = change_edge_num;
    res["diff_edge_num"] = diff_num;
*/
}   

