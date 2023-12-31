#include <iostream>
#include <queue>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <omp.h>
#include <bits/stdc++.h>
#include <iostream>
#include <sys/stat.h>
#include <sys/types.h>
#include <random>

#include "Network/Graph/Graph.h"
#include "Algorithm/AlgorithmBase/AlgorithmBase.h"
#include "Algorithm/Greedy/Greedy.h"
#include "Algorithm/QCAST/QCAST.h"
#include "Algorithm/REPS/REPS.h"
// #include "Algorithm/MyAlgo3/MyAlgo3.h"
// #include "Algorithm/MyGreedyAlgo/MyGreedyAlgo.h"
#include "Algorithm/MyAlgo5/MyAlgo5.h"
using namespace std;




Request generate_new_request(int num_of_node, int time_limit, int min_request, int max_request, int min_value, int max_value){
    //亂數引擎 
    random_device rd;
    default_random_engine generator = default_random_engine(rd());
    uniform_int_distribution<int> unif1(0, num_of_node-1);
    int node1 = unif1(generator), node2 = unif1(generator);
    while(node1 == node2) node2 = unif1(generator);
    
    uniform_int_distribution<int> unif2(min_request, max_request);
    int request = unif2(generator);

    uniform_int_distribution<int> unif3(min_value, max_value);
    int value = unif3(generator);
    return Request(node1, node2, time_limit, request, value);
}
/*
Request generate_fix_request(int node1, int node2, int time_limit, int request){//demo
    return Request(node1, node2, time_limit, request);
}
*/


void create_dir_if_not_exists(const std::string &path) {
	struct stat info;
	if (stat(path.c_str(), &info) == 0 && info.st_mode & S_IFDIR) {
		return;
	}
	mkdir(path.c_str(), 0777);
	return;
}


int main(int argc, char *argv[]){
    string file_path = "../data/";
    if(argc > 1){
	cout<<argv[1]<<endl;
	file_path = string(argv[1]) + '/';
    	cerr<<"the data is store in "<<file_path<<endl;
	create_dir_if_not_exists(file_path);
	create_dir_if_not_exists(file_path+"ans");
	create_dir_if_not_exists(file_path+"input");
	create_dir_if_not_exists(file_path+"log");
    }
    map<string, double> default_setting;
    default_setting["num_of_node"] = 70;
    default_setting["area_alpha"] = 0.4;
    default_setting["memory_cnt_avg"] = 15;
    default_setting["channel_cnt_avg"] = 7;
    default_setting["new_request_cnt"] = 40;
    default_setting["resource_ratio"] = 1;
    default_setting["swap_prob"] = 0.9;
    default_setting["entangle_alpha"] = 0.0002;
    default_setting["total_time_slot"] = 1;
    default_setting["request_avg"] = 3;
    default_setting["epsilon"] = 0.2;    
    default_setting["value"] = 30;
    default_setting["given_path_num"] = 10 ;
    // not used in this paper
    default_setting["node_time_limit"] = 1;
    default_setting["social_density"] = 0.5;
    default_setting["min_fidelity"] = 0.7;
    default_setting["max_fidelity"] = 0.95;
    default_setting["swap_ratio"] = 1; 
    default_setting["request_time_limit"] = 1;
    default_setting["service_time"] = 100;
    map<string, vector<double>> change_parameter;
    change_parameter["swap_prob"] = {0.75, 0.8 , 0.85, 0.9 ,0.95};
    change_parameter["entangle_alpha"] = {0.0004, 0.0003,0.0002, 0.0001, 0};
    change_parameter["min_fidelity"] = {0.5, 0.7, 0.75, 0.85, 0.95};
    change_parameter["resource_ratio"] = {0.5, 1, 1.5, 2, 2.5};
    change_parameter["area_alpha"] = {0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9}; 
    change_parameter["social_density"] = {0.25, 0.5, 0.75, 1}; 
    change_parameter["swap_ratio"] = {0.25, 0.5, 1,  1.5, 2};
    change_parameter["new_request_cnt"] = {20, 30, 40, 50, 60};
    //change_parameter["request_avg"] = {3,4,5};
    change_parameter["num_of_node"] = { /*50,60,70,80,90,*/100};
    change_parameter["memory_cnt_avg"] = { 5 , 7, 9, 11 , 13};
    change_parameter["given_path_num"] = { 10, 11, 12, 13, 14, 15};

    vector<string> X_names =  { /*"new_request_cnt" ,"swap_prob","entangle_alpha","resource_ratio",   "memory_cnt_avg" , "area_alpha" ,"given_path_num",*/"num_of_node"/*,,"request_avg"*/}; 
    vector<string> Y_names =  {  "runtime", "use_channel","drop_req_no" ,"total_earn", "value_per_memory", "value_per_channel","use_memory"
                             /*"throughputs" ,"max_over_ratio","S_D_complete_ratio_difference", "path_success_avg" ,
                             "path_success_avg_before_ent", "new_success_ratio","use_memory_ratio","use_channel_ratio",
			                 "divide_cnt", "change_edge_num", "diff_edge_num", "diff_rate","edge_difference"*/};
    vector<string> algo_names = { /*"MyAlgo3",*/"Greedy_Nonlimit","QCAST_Nonlimit","REPS_Nonlimit",/* "MyAlgo3_0.100000", "MyAlgo3_0.300000,*/ "MyAlgo5"};//{ "MyAlgo3_0.400000","MyAlgo3_0.600000", "MyAlgo3_0.800000"}; //"MyAlgo", "MyGreedyAlgo", "MyAlgo2", 

    // init result
    for(string X_name : X_names) {
        for(string Y_name : Y_names){
            string filename = "ans/" + X_name + "_" + Y_name + ".ans";
            fstream file( file_path + filename, ios::out );
        }
    }
    random_device rd2;
    default_random_engine generator2 = default_random_engine(rd2());;
    std::normal_distribution<double> normal_distribution(1.0,1.5);

    int round = 150;
    for(string X_name : X_names) {
        map<string, double> input_parameter = default_setting;
        for(double change_value : change_parameter[X_name]) {         
            vector<map<string, map<string, double>>> result(round);

            map<string,map<string,vector<double>>> sum_vt;

            input_parameter[X_name] = change_value;
            int num_of_node = input_parameter["num_of_node"];
            int given_path_num = input_parameter["given_path_num"];
            // double social_density = input_parameter["social_density"];
            double area_alpha = input_parameter["area_alpha"];
            double resource_ratio = input_parameter["resource_ratio"];
            int min_memory_cnt = input_parameter["memory_cnt_avg"] * resource_ratio - 2;
            int max_memory_cnt = input_parameter["memory_cnt_avg"] * resource_ratio + 2;
            int min_channel_cnt = input_parameter["channel_cnt_avg"] * resource_ratio - 2;
            int max_channel_cnt = input_parameter["channel_cnt_avg"] * resource_ratio + 2;
            int max_request = input_parameter["request_avg"] + 1;
            int min_request = input_parameter["request_avg"] - 1;
            int max_value = input_parameter["value"] + 5;
            int min_value = input_parameter["value"] - 5;

            double min_fidelity = input_parameter["min_fidelity"];
            double max_fidelity = input_parameter["max_fidelity"];
            double swap_ratio = input_parameter["swap_ratio"];
            double swap_prob = input_parameter["swap_prob"], entangle_alpha = input_parameter["entangle_alpha"];
            double min_swap_prob = input_parameter["swap_prob"] - 0.1;
            double max_swap_prob = input_parameter["swap_prob"] + 0.1;
            int node_time_limit = input_parameter["node_time_limit"];
            int new_request_cnt = input_parameter["new_request_cnt"];
            int service_time = input_parameter["service_time"];
            int request_time_limit = input_parameter["request_time_limit"];
            int total_time_slot = input_parameter["total_time_slot"];
            // python generate graph


            #pragma omp parallel for
            for(int T = 0; T < round; T++){
                string round_str = to_string(T);
                ofstream ofs;
                ofs.open(file_path + "log/" + X_name + "_in_" + to_string(change_value) + "_Round_" + round_str + ".log");

                time_t now = time(0);
                char* dt = ctime(&now);
                cerr  << "時間 " << dt << endl << endl; 
                ofs  << "時間 " << dt << endl << endl; 

                string filename = file_path + "input/round_" + round_str + ".input";
                string command = "python3 main.py ";
                string parameter = to_string(num_of_node) + " " + to_string(min_channel_cnt) + " " + to_string(max_channel_cnt) + " " + to_string(min_memory_cnt) + " " + to_string(max_memory_cnt) + " " + to_string(min_fidelity) + " " + to_string(max_fidelity) + " " + to_string(area_alpha) + " " + to_string(min_swap_prob) + " " +  to_string(max_swap_prob);
                //cout<<command + filename + " " + parameter<<endl;
                if(system((command + filename + " " + parameter).c_str()) != 0){
                    cerr<<"error:\tsystem proccess python error"<<endl;
                    exit(1);
                }

                
                vector<AlgorithmBase*> algorithms;
                //algorithms.emplace_back(new MyAlgo3(filename, request_time_limit, node_time_limit, swap_prob, entangle_alpha));
                //algorithms.emplace_back(new Greedy(filename, request_time_limit, node_time_limit, swap_prob, entangle_alpha , true));
                algorithms.emplace_back(new Greedy(filename, request_time_limit, node_time_limit, swap_prob, entangle_alpha , false));
                //algorithms.emplace_back(new QCAST(filename, request_time_limit, node_time_limit, swap_prob, entangle_alpha , true));
                algorithms.emplace_back(new QCAST(filename, request_time_limit, node_time_limit, swap_prob, entangle_alpha , false));
                //algorithms.emplace_back(new REPS(filename, request_time_limit, node_time_limit, swap_prob, entangle_alpha ,true ));
                algorithms.emplace_back(new REPS(filename, request_time_limit, node_time_limit, swap_prob, entangle_alpha ,false ));
                //algorithms.emplace_back(new MyAlgo3(filename, request_time_limit, node_time_limit, swap_prob, entangle_alpha));
                //algorithms.emplace_back(new MyAlgo3(filename, request_time_limit, node_time_limit, swap_prob, entangle_alpha, 0.1 ));
                //algorithms.emplace_back(new MyAlgo3(filename, request_time_limit, node_time_limit, swap_prob, entangle_alpha, 0.3));
                //algorithms.emplace_back(new MyAlgo3(filename, request_time_limit, node_time_limit, swap_prob, entangle_alpha, 0.4 ));
                //algorithms.emplace_back(new MyAlgo3(filename, request_time_limit, node_time_limit, swap_prob, entangle_alpha, 0.6 ));
                //algorithms.emplace_back(new MyAlgo3(filename, request_time_limit, node_time_limit, swap_prob, entangle_alpha, 0.8 )); 
                algorithms.emplace_back(new MyAlgo5(filename, request_time_limit, node_time_limit, swap_prob, entangle_alpha , given_path_num));
                // 建完圖，刪除 input 檔避免佔太多空間
                // command = "rm -f " + file_path + "input/round_" + round_str + ".input";
                // if(system((command).c_str()) != 0){
                //     cerr<<"error:\tsystem proccess delete input file error"<<endl;
                //     exit(1);
                // }
                ofs<<"---------------in round " <<T<<" -------------" <<endl;
                for(int t = 0; t < total_time_slot; t++){
                    ofs<<"---------------in timeslot " <<t<<" -------------" <<endl;
                    cout<< "---------generating requests in main.cpp----------" << endl;

                    for(int q = 0; q < new_request_cnt && t < service_time; q++){
                        bool check_no_repeat;
                        do{
                            check_no_repeat=true;
                            Request new_request = generate_new_request(num_of_node, request_time_limit, min_request, max_request, min_value, max_value);
                            //cout<<"new"<<new_request.get_source()<<"----------->"<<new_request.get_destination()<<endl;

                            for(auto &req_num:algorithms[0]->get_requests()){
                                //cout<<req_num.get_source()<<"------------->"<<req_num.get_destination()<<endl;
                                if(req_num.get_source()==new_request.get_source() && req_num.get_destination()==new_request.get_destination()){
                                    check_no_repeat=false;
                                }
                            }
                            if(check_no_repeat==true){
                                cout<<q << ". source: " << new_request.get_source()<<", destination: "<<new_request.get_destination()<<endl;

                                //------------------setting value to requests------------------
                                //algorithms[0]->add_new_request(new_request);
                            
                            
                                double willness=0;
                                while(willness<1.0 || willness>3.0){
                                    willness = normal_distribution(generator2);
                                }
                                new_request.set_willness(willness);
                                //cout<<"willness:"<<new_request.get_willness()<<endl;

                                //--------------------graph------------------------------------
                                ifstream graph_input;
                                int num_of_node_copy;
                                int num_of_edge_copy;
                                graph_input.open (filename);
                                graph_input >> num_of_node_copy;
                                vector<vector<int>> neighbor;
                                neighbor.resize(num_of_node_copy);
                                vector<int>(num_of_node_copy,-1);

                                double pos_x, pos_y, new_swap_prob;
                                int memory_cnt;
                                for(int i = 0; i < num_of_node; i++){
                                    graph_input >> pos_x >> pos_y >> memory_cnt >> new_swap_prob;
                                }
                                int node_id1, node_id2;
                                int channel_cnt;
                                double fidelity;
                                double dis_sum = 0;
                                double prob_sum = 0;
                                graph_input >> num_of_edge_copy;
                                for(int i = 0;i < num_of_edge_copy; i++){
                                    graph_input >> node_id1 >> node_id2 >> channel_cnt >> fidelity;
                                    neighbor[node_id1].emplace_back(node_id2);
                                    neighbor[node_id2].emplace_back(node_id1);
                                }
                                //-----------------BFS-------------------------
                                vector<int>color(num_of_node_copy,0);
                                //vector<int>predecessor(num_of_node_copy,-1);
                                vector<int>distance(num_of_node_copy,num_of_node_copy+1);

                                std::queue<int> queue;
                                int src = new_request.get_source();

                                for (int j = 0; j < num_of_node_copy; j++) {  // j從0數到num_of_node_copy-1, 因此j會走過graph中所有vertex
                                    if (color[src] == 0) {                // 第一次i會是起點vertex, 如圖二(c)
                                        color[src] = 1;                   // 1:灰色
                                        distance[src] = 0;                // 每一個connected component的起點之距離設成0
                                        //predecessor[i] = -1;            // 每一個connected component的起點沒有predecessor
                                        queue.push(src);
                                        while (!queue.empty()) {
                                            int u = queue.front();                  // u 為新的搜尋起點
                                            for (auto neigh:neighbor[u]) {                         // 分成兩段
                                                if (color[neigh] == 0) {                // 若被「找到」的vertex是白色
                                                    color[neigh] = 1;                   // 塗成灰色, 表示已經被「找到」
                                                    distance[neigh] = distance[u] + 1;  // 距離是predecessor之距離加一
                                                    //predecessor[neigh] = u;             // 更新被「找到」的vertex的predecessor
                                                    queue.push(neigh);                      // 把vertex推進queue
                                                }
                                            }
                                            queue.pop();        // 把u移出queue
                                            color[u] = 2;   // 並且把u塗成黑色
                                        }
                                    }
                                    // 若一次回圈沒有把所有vertex走過, 表示graph有多個connected component
                                    // 就把i另成j, 繼續檢查graph中的其他vertex是否仍是白色, 若是, 重複while loop
                                    if(distance[new_request.get_destination()]!=num_of_node_copy+1){
                                        break;
                                    }
                                    src = j;
                                }
                                new_request.set_value(distance[new_request.get_destination()]+distance[new_request.get_destination()]*swap_ratio);
                                //cout<<"new_request "<<q<<" with value:"<<new_request.get_value()<<" and willness:"<<new_request.get_willness()<<" in round "<<T<<endl;
                                for(auto &algo:algorithms){
                                    result[T][algo->get_name()]["total_request"]++; 
                                    algo->add_new_request(new_request);
                                    //cout<<"what happend "<<algorithms[0]->get_requests().size()<<endl;
                                }
                                graph_input.close();
                            }   
                        }while(check_no_repeat==false);
                    }

                    cout<< "---------generating requests in main.cpp----------end" << endl;
                    

                    //#pragma omp parallel for 
                    for(int i = 0; i < (int)algorithms.size(); i++){
                        auto &algo = algorithms[i];
                        ofs<<"-----------run "<< algo->get_name() << " ---------"<<endl;
                        algo->run();

                        //ofs<<"total_throughputs : "<<algo->get_res("throughputs")<<endl;
                        ofs<<"total_earn :"<<algo->get_res("total_earn")<<endl;
                        ofs<<"-----------run "<<algo->get_name() << " ---------end"<<endl;
                    }
                    
                }
                ofs<<"---------------in round " <<T<<" -------------end" <<endl;
                ofs << endl;
                for(auto &algo:algorithms){
                    //ofs<<"("<<algo->get_name()<<")total throughput = "<<algo->get_res("throughputs")<<endl;
                    ofs<<"("<<algo->get_name()<<")total earn = "<<algo->get_res("total_earn")<<endl;
                }
                cout<<"---------------in round " <<T<<" -------------end" <<endl;
                cout << endl;
                for(auto &algo:algorithms){
                    //cout<<"("<<algo->get_name()<<")total throughput = "<<algo->get_res("throughputs")<<endl;
                    cout<<fixed<<"("<<algo->get_name()<<")total earn = "<<algo->get_res("total_earn")<<endl;
                }
                
                for(auto &algo:algorithms){
                    if(algo->get_name()=="MyAlgo5"){
                        result[T]["MyAlgo5"]["UB"]=algo->get_res("UB");
                    }
                    for(string Y_name : Y_names) {
                        result[T][algo->get_name()][Y_name] = algo->get_res(Y_name);
                    }
                }

                

                
                for(auto &algo:algorithms){
                    for(string Y_name :Y_names){
                        for(auto it:algo->get_res_vt()){
                            sum_vt[algo->get_name()][Y_name].push_back(it);
                        }
                    }
                }


                now = time(0);
                dt = ctime(&now);
                cerr  << "時間 " << dt << endl << endl; 
                ofs  << "時間 " << dt << endl << endl; 
                ofs.close();
            
                for(auto &algo:algorithms){
                    delete algo;
                }
                algorithms.clear();
            
            }
            
            map<string, map<string, double>> sum_res;
             for(string algo_name : algo_names){
                 for(int T = 0; T < round; T++){
            //         result[T][algo_name]["waiting_time"] /= result[T][algo_name]["total_request"];
            //         result[T][algo_name]["encode_ratio"] = result[T][algo_name]["encode_cnt"] / (result[T][algo_name]["encode_cnt"] + result[T][algo_name]["unencode_cnt"]);
            //         result[T][algo_name]["succ-finished_ratio"] = result[T][algo_name]["throughputs"] / result[T][algo_name]["finished_throughputs"];
            //         result[T][algo_name]["fail-finished_ratio"] = 1 - result[T][algo_name]["succ-finished_ratio"];
            //         result[T][algo_name]["path_length"] = result[T][algo_name]["path_length"] / result[T][algo_name]["finished_throughputs"];
            //         result[T][algo_name]["divide_cnt"] = result[T][algo_name]["divide_cnt"] / result[T][algo_name]["finished_throughputs"];
            //        result[T][algo_name]["use_memory_ratio"] = result[T][algo_name]["use_memory"] / result[T][algo_name]["total_memory"];
            //       result[T][algo_name]["use_channel_ratio"] = result[T][algo_name]["use_channel"] / result[T][algo_name]["total_channel"];
                     result[T][algo_name]["value_per_memory"] = result[T][algo_name]["total_earn"] / result[T][algo_name]["use_memory"]  ;                  
                     result[T][algo_name]["value_per_channel"] = result[T][algo_name]["total_earn"] / result[T][algo_name]["use_channel"];
                 }
             }

            /*
            for(int T = 0; T < round; T++){
                // result[T]["MyAlgo3"]["diff_rate"] = result[T]["MyAlgo3"]["change_edge_num"] / result[T]["MyAlgo3"]["diff_edge_num"];
                result[T]["MyAlgo"]["edge_difference"] = result[T]["MyAlgo"]["change_edge_num"] - result[T]["MyAlgo3"]["change_edge_num"];
            }
            */

            for(int T = 0; T < round; T++){
            //cout<<result[T]["MyAlgo3"]["primal"]<<" "<<result[T]["MyAlgo3_0.100000"]["primal"]<<" "<<result[T]["MyAlgo3_0.300000"]["primal"]<<endl;
              sum_res["MyAlgo5"]["UB"] += result[T]["MyAlgo5"]["UB"];
              //cout<<"Outside UB:"<< result[T]["MyAlgo5"]["UB"]<<endl;
            }
            
                

            for(string Y_name : Y_names) {
                string filename = "ans/" + X_name + "_" + Y_name + ".ans";
                ofstream ofs;
                ofs.open(file_path + filename, ios::app);
                ofs << change_value << ' ';
                
                for(string algo_name : algo_names){
                    for(int T = 0; T < round; T++){
                        sum_res[algo_name][Y_name] += result[T][algo_name][Y_name];
                    }
                    ofs << sum_res[algo_name][Y_name] / round << ' ';
                    
                }
                
                if(Y_name == "total_earn"){
                    ofs << sum_res["MyAlgo5"]["UB"] / round << " ";
                }
                
                ofs << endl;
                ofs.close();
            }

            

                
            // string filename = "ans/" + X_name + "_success_path_prob_vt.ans";
            // ofstream ofs;
            // ofs.open(file_path + filename, ios::app);
            // ofs << change_value << endl;
            
            // for(string algo_name : algo_names){
            //     ofs<<algo_name<<endl;
            //     sort(sum_vt[algo_name]["runtime"].begin(),sum_vt[algo_name]["runtime"].end());
            //     for(auto it:sum_vt[algo_name]["runtime"]){
            //         ofs << it << " ";
            //     }
            //     ofs << endl;
            // }
            
            // ofs << endl;
            // ofs.close();
            
            
        }
    }
    return 0;
}
