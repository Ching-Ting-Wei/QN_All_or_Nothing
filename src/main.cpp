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
    default_setting["area_alpha"] = 0.6;
    default_setting["memory_cnt_avg"] = 12;
    default_setting["channel_cnt_avg"] = 5;
    default_setting["new_request_cnt"] = 30;
    default_setting["resource_ratio"] = 1;
    default_setting["swap_prob"] = 0.97;
    default_setting["entangle_alpha"] = 0.0002;
    default_setting["total_time_slot"] = 1;
    default_setting["request_avg"] = 3;
    default_setting["epsilon"] = 0.1;    
    default_setting["value"] = 25;
    default_setting["given_path_num"] = 10;
    // not used in this paper
    default_setting["node_time_limit"] = 1;
    default_setting["social_density"] = 0.5;
    default_setting["min_fidelity"] = 0.7;
    default_setting["max_fidelity"] = 0.95;
    default_setting["request_time_limit"] = 1;
    default_setting["service_time"] = 100;

    map<string, vector<double>> change_parameter;
    change_parameter["swap_prob"] = {0.75, 0.8 , 0.85, 0.9 ,0.95};
    change_parameter["entangle_alpha"] = {0.01 ,0.005 ,0.001, 0.0005, 0.0001};
    change_parameter["min_fidelity"] = {0.5, 0.7, 0.75, 0.85, 0.95};
    change_parameter["resource_ratio"] = {0.5, 1, 1.5, 2, 2.5};
    change_parameter["area_alpha"] = {0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9}; 
    change_parameter["social_density"] = {0.25, 0.5, 0.75, 1}; 
    change_parameter["new_request_cnt"] = {10,20,30,40,50};
    //change_parameter["request_avg"] = {3,4,5};
    change_parameter["num_of_node"] = { 50, 60, 70 , 80, 90, 100};
    change_parameter["memory_cnt_avg"] = { 5 , 7, 9, 11 , 13};
    change_parameter["given_path_num"] = { 10, 11, 12, 13, 14, 15};

    vector<string> X_names =  { "entangle_alpha" }; 
    vector<string> Y_names =  {  "use_memory", "use_channel","drop_req_no" ,"total_earn", "value_per_memory", "value_per_channel"};
    vector<string> algo_names = { /*"MyAlgo3",*/"Greedy_Nonlimit","QCAST_Nonlimit","REPS_Nonlimit",/* "MyAlgo3_0.100000", "MyAlgo3_0.300000,*/ "MyAlgo5"};//{ "MyAlgo3_0.400000","MyAlgo3_0.600000", "MyAlgo3_0.800000"}; //"MyAlgo", "MyGreedyAlgo", "MyAlgo2", 

    // init result

    for(string Y_name : Y_names){
        string filename = "ans/entangle_alpha_" + Y_name + ".ans";
        fstream file( file_path + filename, ios::out );
    }

    int round = 50;
    map<string, double> input_parameter = default_setting;
    vector<vector<Request>>request_list(round);
    vector<map<string, map<string, double>>> result(round);
    vector<vector<map<string, map<string, double>>>> result2(round);
    for(auto &it:result2){
        it.resize(20);
    }
    map<string,map<string,vector<double>>> sum_vt;
    
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
    int max_value = input_parameter["value"] + 10;
    int min_value = input_parameter["value"] - 10;
    double min_fidelity = input_parameter["min_fidelity"];
    double max_fidelity = input_parameter["max_fidelity"];

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
    for(int T = 0; T < round; T++){           //做T次圖 
        string round_str = to_string(T);
        string filename = file_path + "input/round_" + round_str + ".input";
        string command = "python3 main.py ";
        string parameter = to_string(num_of_node) + " " + to_string(min_channel_cnt) + " " + to_string(max_channel_cnt) + " " + to_string(min_memory_cnt) + " " + to_string(max_memory_cnt) + " " + to_string(min_fidelity) + " " + to_string(max_fidelity) + " " + to_string(area_alpha) + " " + to_string(min_swap_prob) + " " +  to_string(max_swap_prob);
        //cout<<command + filename + " " + parameter<<endl;
        if(system((command + filename + " " + parameter).c_str()) != 0){
            cerr<<"error:\tsystem proccess python error"<<endl;
            exit(1);
        }
        for(int q = 0; q < new_request_cnt; q++){
            bool check_no_repeat;
            do{
                check_no_repeat=true;
                Request new_request = generate_new_request(num_of_node, request_time_limit, min_request, max_request, min_value, max_value);
                for(auto &it:request_list[T]){
                    if(it.get_source()==new_request.get_source() && it.get_destination()==new_request.get_destination()){
                        check_no_repeat=false;
                    }
                }
                if(check_no_repeat==true){
                    request_list[T].push_back(new_request);
                }
            }while(check_no_repeat==false);
        }
        
    }


    //#pragma omp parallel for
    for(int change_index = 0;change_index < change_parameter["entangle_alpha"].size();change_index++) { 
        ofstream ofs;
        time_t now;
        char* dt;
        double no = change_parameter["entangle_alpha"][change_index];
        now = time(0);;
        dt = ctime(&now);
        #pragma omp parallel for
        for(int T = 0; T < round; T++){
            string round_str = to_string(T);
            cout<<"entangle_alpha:"<<no<<" round:"<<round_str<<"-------------------------"<<endl;
            ofs.open(file_path + "log/" + "entangle_alpha" + "_in_" + to_string(no) + "_Round_" + round_str + ".log");
            cerr  << "時間 " << dt << endl << endl; 
            ofs  << "時間 " << dt << endl << endl; 
            
            
            string filename = file_path + "input/round_" + round_str + ".input";
            vector<AlgorithmBase*> algorithms;
            algorithms.emplace_back(new Greedy(filename, request_time_limit, node_time_limit, swap_prob, no , false));
            algorithms.emplace_back(new QCAST(filename, request_time_limit, node_time_limit, swap_prob, no , false));
            algorithms.emplace_back(new REPS(filename, request_time_limit, node_time_limit, swap_prob, no ,false ));
            algorithms.emplace_back(new MyAlgo5(filename, request_time_limit, node_time_limit, swap_prob, no , given_path_num));

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
                    int q=0;
                    for(auto &req:request_list[T]){
                        if(t >= service_time){break;}
                        Request new_request = req;
                        cout<< q << ". source: " << new_request.get_source()<<", destination: "<<new_request.get_destination()<<endl;
                        for(auto &algo:algorithms){
                            result[T][algo->get_name()]["total_request"]++; 
                            algo->add_new_request(new_request);
                        }
                        q++;
                    }

                    cout<< "---------generating requests in main.cpp----------end" << endl;
                    #pragma omp parallel for 
                    for(int i = 0; i < (int)algorithms.size(); i++){
                        auto &algo = algorithms[i];
                        ofs<<"-----------run "<< algo->get_name() << " ---------"<<endl;
                        algo->run();
                        ofs<<"total_earn :"<<algo->get_res("total_earn")<<endl;
                        ofs<<"-----------run "<<algo->get_name() << " ---------end"<<endl;
                    }
                }
                
                ofs<<"---------------in round " <<T<<" -------------end" <<endl;
                ofs << endl;
                cout<<"---------------in round " <<T<<" -------------end" <<endl;
                cout << endl;

                for(auto &algo:algorithms){
                    ofs<<"("<<algo->get_name()<<")total earn = "<<algo->get_res("total_earn")<<endl;
                }

                for(auto &algo:algorithms){
                    //cout<<"("<<algo->get_name()<<")total throughput = "<<algo->get_res("throughputs")<<endl;
                    cout<<fixed<<"("<<algo->get_name()<<")total earn = "<<algo->get_res("total_earn")<<endl;
                }
                
                for(auto &algo:algorithms){
                    for(string Y_name : Y_names) {
                        result[T][algo->get_name()][Y_name] = algo->get_res(Y_name);
                    }
                }
                cout<<"there>>>>>>>>>>>>>>>>>>>"<<endl;

                for(auto &algo:algorithms){  
                    for(string Y_name : Y_names) {  
                        result2[T][change_index][algo->get_name()][Y_name] = algo->get_res(Y_name);
                    }
                    //cout<<"index:"<<change_index<<" "<<T<<endl;
                    //cout<<"value:"<<result2[T][change_index][algo->get_name()]["new_request_cnt"]<<endl;
                }

                cout<<"there2>>>>>>>>>>>>>>>>>>>"<<endl;
                now = time(0);;
                dt = ctime(&now);
                cerr  << "時間 " << dt << endl << endl; 
                ofs  << "時間 " << dt << endl << endl; 
                ofs.close();
            
                for(auto &algo:algorithms){
                    delete algo;
                }
                algorithms.clear();
        }
        cout<<"here>>>>>>>>>>>>>>>>>>>>>>>>>>"<<endl;
        vector<map<string, map<string, double>>> sum_res(10);

        for(string algo_name : algo_names){
            for(int T = 0; T < round; T++){
                for(int index = 0;index <change_parameter["entangle_alpha"].size();index++){
                    result2[T][change_index][algo_name]["value_per_memory"] = result2[T][change_index][algo_name]["total_earn"] / result2[T][change_index][algo_name]["use_memory"]  ;                  
                    result2[T][change_index][algo_name]["value_per_channel"] = result2[T][change_index][algo_name]["total_earn"] / result2[T][change_index][algo_name]["use_channel"];
                }

            }
        }
        
        for(string Y_name : Y_names) {
            string filename = "ans/entangle_alpha_" + Y_name + ".ans";
            ofstream ofs;
            ofs.open(file_path + filename, ios::app);
            ofs << change_parameter["entangle_alpha"][change_index] << ' ';
        
            for(string algo_name : algo_names){

                    for(int T = 0; T < round; T++){
                        sum_res[change_index][algo_name][Y_name] += result2[T][change_index][algo_name][Y_name];
                        cout<<change_index<<"------value:"<<result2[T][change_index][algo_name][Y_name]<<endl;
                    }
                    ofs << sum_res[change_index][algo_name][Y_name] / round << ' ';
                
            }
            /*
            if(Y_name == "throughputs"){
                ofs << sum_res["MyAlgo3"]["primal"] / round << " ";
            }
            */
            ofs << endl;
            ofs.close();
        }
    }

    return 0;
}
