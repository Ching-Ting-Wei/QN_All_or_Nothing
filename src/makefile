gurobi_path = /home/oba/danger/gurobi951/linux64/
g++ = g++ -std=c++17 -fopenmp -g -Wall -Wextra -O3

# on server
#gurobi_path=/home/oba/danger/gurobi951/linux64/

dependency = Network/Node/Node.o Network/Channel/Channel.o Path/Path.o Network/Graph/Graph.o Algorithm/AlgorithmBase/AlgorithmBase.o Request/Request.o config.o

all: a.out
encode_path.out: encode_path_ratio.cpp Network/Node/Node.o Network/Channel/Channel.o Path/Path.o Network/Graph/Graph.o Algorithm/AlgorithmBase/AlgorithmBase.o Request/Request.o config.o Algorithm/MyAlgo/MyAlgo.o Request/WholeRequest.o
	$(g++) encode_path_ratio.cpp Network/Node/Node.o Network/Channel/Channel.o Path/Path.o Network/Graph/Graph.o Algorithm/AlgorithmBase/AlgorithmBase.o Request/Request.o config.o Algorithm/MyAlgo/MyAlgo.o Request/WholeRequest.o -o encode_path.out

a.out: main.cpp $(dependency) Algorithm/Greedy/Greedy.o Request/Request.o Algorithm/MyAlgo/MyAlgo.o Algorithm/MyAlgo2/MyAlgo2.o  Algorithm/MyAlgo3/MyAlgo3.o Algorithm/MyAlgo5/MyAlgo5.o Algorithm/QCAST/QCAST.o Algorithm/REPS/REPS.o Algorithm/MyAlgo5/PathGraph.o
	$(g++) main.cpp $(dependency) Algorithm/Greedy/Greedy.o Algorithm/MyAlgo/MyAlgo.o Algorithm/MyAlgo2/MyAlgo2.o Algorithm/MyAlgo3/MyAlgo3.o Algorithm/MyAlgo5/MyAlgo5.o Algorithm/QCAST/QCAST.o Algorithm/REPS/REPS.o Algorithm/MyAlgo5/PathGraph.o -I$(gurobi_path)/include/ -L$(gurobi_path)/lib -lgurobi_c++ -lgurobi95 -lm 

#fight_with_local_opt.out: fight_with_local_opt.cpp $(dependency) Request/WholeRequest.o 
#	$(g++) fight_with_local_opt.cpp $(dependency) Algorithm/MyAlgo/MyAlgo.o  Algorithm/MyAlgo2/MyAlgo2.o Algorithm/MyGreedyAlgo/MyGreedyAlgo.o Request/WholeRequest.o -o fight_with_local_opt.out

# qcast.out: Network/Node/Node.o Network/Channel/Channel.o Path/Path.o Network/Graph/Graph.o Algorithm/AlgorithmBase/AlgorithmBase.o
Network/Node/Node.o: Network/Node/Node.cpp Network/Node/Node.h config.h
	$(g++) -c Network/Node/Node.cpp Network/Node/Node.h config.h
	mv Node.o Network/Node/Node.o

Network/Channel/Channel.o: Network/Channel/Channel.h Network/Channel/Channel.cpp config.h
	$(g++) -c Network/Channel/Channel.h Network/Channel/Channel.cpp config.h
	mv Channel.o Network/Channel/Channel.o 

Path/Path.o: Path/Path.cpp Path/Path.h config.h
	$(g++) -c Path/Path.cpp Path/Path.h config.h
	mv Path.o Path/Path.o

Network/Graph/Graph.o: Network/Graph/Graph.cpp Network/Graph/Graph.h config.h
	$(g++) -c Network/Graph/Graph.cpp Network/Graph/Graph.h config.h
	mv Graph.o Network/Graph/Graph.o

Request/Request.o: Request/Request.cpp Request/Request.h config.h
	$(g++) -c Request/Request.cpp Request/Request.h config.h
	mv Request.o Request/Request.o
	
Request/WholeRequest.o: Request/WholeRequest.cpp Request/WholeRequest.h config.h
	$(g++) -c Request/WholeRequest.cpp Request/WholeRequest.h config.h
	mv WholeRequest.o Request/WholeRequest.o

Algorithm/AlgorithmBase/AlgorithmBase.o: Algorithm/AlgorithmBase/AlgorithmBase.h Algorithm/AlgorithmBase/AlgorithmBase.cpp config.h
	$(g++) -c Algorithm/AlgorithmBase/AlgorithmBase.h Algorithm/AlgorithmBase/AlgorithmBase.cpp config.h
	mv AlgorithmBase.o Algorithm/AlgorithmBase/AlgorithmBase.o

Algorithm/Greedy/Greedy.o: Algorithm/Greedy/Greedy.h Algorithm/Greedy/Greedy.cpp config.h
	$(g++) -c Algorithm/Greedy/Greedy.h Algorithm/Greedy/Greedy.cpp config.h
	mv Greedy.o Algorithm/Greedy/Greedy.o

Algorithm/MyAlgo/MyAlgo.o: Algorithm/MyAlgo/MyAlgo.h Algorithm/MyAlgo/MyAlgo.cpp config.h
	$(g++) -c --static Algorithm/MyAlgo/MyAlgo.h Algorithm/MyAlgo/MyAlgo.cpp config.h
	mv MyAlgo.o Algorithm/MyAlgo/MyAlgo.o

Algorithm/MyAlgo2/MyAlgo2.o: Algorithm/MyAlgo2/MyAlgo2.h Algorithm/MyAlgo2/MyAlgo2.cpp config.h
	$(g++) -c --static Algorithm/MyAlgo2/MyAlgo2.h Algorithm/MyAlgo2/MyAlgo2.cpp config.h
	mv MyAlgo2.o Algorithm/MyAlgo2/MyAlgo2.o

Algorithm/MyAlgo3/MyAlgo3.o: Algorithm/MyAlgo3/MyAlgo3.h Algorithm/MyAlgo3/MyAlgo3.cpp config.h
	$(g++) -c --static Algorithm/MyAlgo3/MyAlgo3.h Algorithm/MyAlgo3/MyAlgo3.cpp config.h
	mv MyAlgo3.o Algorithm/MyAlgo3/MyAlgo3.o

Algorithm/MyAlgo5/MyAlgo5.o: Algorithm/MyAlgo5/MyAlgo5.h Algorithm/MyAlgo5/MyAlgo5.cpp config.h
	$(g++) -c --static Algorithm/MyAlgo5/MyAlgo5.h Algorithm/MyAlgo5/MyAlgo5.cpp config.h
	mv MyAlgo5.o Algorithm/MyAlgo5/MyAlgo5.o

Algorithm/QCAST/QCAST.o: Algorithm/QCAST/QCAST.h Algorithm/QCAST/QCAST.cpp config.h
	$(g++) -c Algorithm/QCAST/QCAST.h Algorithm/QCAST/QCAST.cpp config.h
	mv QCAST.o Algorithm/QCAST/QCAST.o

Algorithm/REPS/REPS.o: Algorithm/REPS/REPS.h Algorithm/REPS/REPS.cpp config.h
	$(g++) -c --static Algorithm/REPS/REPS.h Algorithm/REPS/REPS.cpp config.h -m64 -I$(gurobi_path)/include/  -L$(gurobi_path)/lib -lgurobi_c++ -lgurobi95 -lm 
	mv REPS.o Algorithm/REPS/REPS.o

Algorithm/MyAlgo5/PathGraph.o: Algorithm/MyAlgo5/PathGraph.h Algorithm/MyAlgo5/PathGraph.cpp config.h
	$(g++) -c --static Algorithm/MyAlgo5/PathGraph.h Algorithm/MyAlgo5/PathGraph.cpp config.h
	mv PathGraph.o Algorithm/MyAlgo5/PathGraph.o

config.o:	config.h config.cpp
	g++ -std=c++17 -c -g -O3 config.h config.cpp
check_node: Node/Node.o Node/Node_test.cpp
	g++ -std=c++17 -g -Wall -Wextra -O3 Node/Node_test.cpp Node/Node.o -o node_check.out

check_channel: Node/Node.o Channel/Channel.o Channel/Channel_test.cpp
	g++ -std=c++17 -g -Wall -Wextra -O3 Node/Node.o Channel/Channel.o Channel/Channel_test.cpp -o channel_check.out


clean:
	rm -f a.out check_node.out channel_check.out
	rm -f Network/Node/Node.o Network/Node/Node.h.gch
	rm -f Network/Channel/Channel.o Network/Channel/Channel.h.gch
	rm -f Path/Path.o Path/Path.h.gch
	rm -f Network/Graph/Graph.o Network/Graph/Graph.h.gch
	rm -f Request/Request.o Request/Request.h.gch
	rm -f Algorithm/MyAlgo/MyAlgo.o Algorithm/MyAlgo/MyAlgo.h.gch
	rm -f Algorithm/MyAlgo2/MyAlgo2.o Algorithm/MyAlgo2/MyAlgo2.h.gch
	rm -f Algorithm/MyAlgo3/MyAlgo3.o Algorithm/MyAlgo3/MyAlgo3.h.gch
	rm -f Algorithm/MyAlgo5/MyAlgo5.o Algorithm/MyAlgo5/MyAlgo5.h.gch
	rm -f Algorithm/MyAlgo5/PathGraph.o Algorithm/MyAlgo5/PathGraph.h.gch
	rm -f config.o config.h.gch
	find . -type f -name '*.o' -delete
	find . -type f -name '*.h.gch' -delete
	# rm -f ../data/ans/*
	# rm -f ../data/log/*
	# rm -f ../data/input/*
