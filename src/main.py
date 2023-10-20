from math import ceil
import sys
import networkx as nx
import random


G = nx.Graph()
vis = {}

def dfs(node):
    global vis
    global G
    if vis[node] == True:
        return
    vis[node] = True
    for v in G.neighbors(node):
        if vis[v] == False:
            dfs(v)
    return

def is_connect():
    global vis
    global G
    #print(type(G.nodes()))
    vis = {node : False for node in G.nodes()}
    for node in G.nodes():
        dfs(node)
        break
    for node in G.nodes():
        if(vis[node] == False):
            return False
    return True

def dist(a, b):
    (p1, p2) = a
    (q1, q2) = b
    return (q1-p1) * (q1-p1) + (q2-p2) * (q2-p2)

def genSocialNetwork(userNum, density):
	# n * n
	community0 = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19] #0
	community1 = [0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 2, 2, 2, 2, 3, 2, 2, 2, 3, 2]  # 0.25
	community2 = [0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1]  # 0.50
	community3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0]  # 0.75
	community4 = [0 for _ in range(20)]                                        # 1.00
	community = {0.0 : community0, 0.25 : community1, 0.50 : community2, 0.75 : community3, 1.00 : community4}
	SN = {i: [] for i in range(userNum)}  # user to user
	for i in range(userNum):
		for j in range(userNum):
			SN[i].append(0)
	community = community[density]
	for i in range(userNum):
		for j in range(i, userNum):
			owni = i % 20
			ownj = j % 20
			if community[owni] == community[ownj]:
				SN[i][j] = 1
				SN[j][i] = 1
			if density == 0.0:
				SN[i][j] = 0
				SN[j][i] = 0
	return SN


if  len(sys.argv) <= 9:
    print("missing argv")
    sys.exit()


filename = sys.argv[1]
num_of_node = int(sys.argv[2])
min_channel_cnt = int(sys.argv[3])
max_channel_cnt = int(sys.argv[4])
min_memory_cnt = int(sys.argv[5])
max_memory_cnt = int(sys.argv[6])
min_fidelity = float(sys.argv[7])
max_fidelity = float(sys.argv[8])
area_alpha = float(sys.argv[9])
min_swap_prob = float(sys.argv[10])
max_swap_prob = float(sys.argv[11])


print("--------------generating graph--------------", file=sys.stderr)
print("--------------generating graph--------------")
print("filename = ", filename, file=sys.stderr)
print("number_of_node = ", num_of_node, file=sys.stderr)
print("min_channel_cnt = ", min_channel_cnt, ", max_channel_cnt = ", max_channel_cnt, file=sys.stderr)
print("min_memory_cnt = ", min_memory_cnt, ", max_memory_cnt = ", max_memory_cnt, file=sys.stderr)
# print("socail density = ", socail_density, file=sys.stderr)


while True:
    G = nx.waxman_graph(num_of_node, beta=0.85, alpha=area_alpha, domain=(0, 0, 0.5, 1))
    #G = nx.waxman_graph(num_of_node, beta=0.85, alpha=10, domain=(0, 0, 0.5, 1))

    positions = nx.get_node_attributes(G, 'pos')
    add_edge = []
    for u in range(G.order()-1):
        mi_dist = dist(positions[u], positions[G.order()-1])
        mi_idx = G.order()-1
        for v in range(u+1, G.order()):
            if(G.has_edge(u, v)):
                continue
            if(mi_dist > dist(positions[u], positions[v])):
                mi_dist = dist(positions[u], positions[v])
                mi_idx = v
        if(G.has_edge(u, mi_idx)):
            continue
        add_edge.append((u, mi_idx))
    for e in add_edge:
        (u, v) = e
        G.add_edge(u, v)
    if is_connect():
        break 
    else:
        print("topo is not connected", file = sys.stderr)

path = filename             #get from argv
f = open(path, 'w') 
positions = nx.get_node_attributes(G, 'pos')
print(num_of_node, file=f)
for n in G.nodes():
    (p1, p2) = positions[n]
    print(str(p1 * 2000) + ' ' + str(p2 * 2000) + " " + str(ceil(random.random()*(max_memory_cnt-min_memory_cnt)) + min_memory_cnt) + " " + str(min(1, random.random()*(max_swap_prob - min_swap_prob) + min_swap_prob)), file=f)
cnt = 0
for e in G.edges():
    if e[0] != e[1]:#file
        cnt += 1
print(cnt, file=f)
for e in G.edges():
    if e[0] != e[1]:#file
        print(str(e[0]) + ' ' + str(e[1]) + " " + str(ceil(random.random()*(max_channel_cnt-min_channel_cnt)) + min_channel_cnt), 1, file=f)


# SN = genSocialNetwork(num_of_node, socail_density)
# for i in range(num_of_node):
#     for j in range(num_of_node):
#         print(SN[i][j], end=" ", file=f)
#     print("", file=f)


#random social
# for i in range(num_of_node):
#     for j in range(num_of_node):
#         if i == j:
#             print(1, end=" ", file=f)
#         else:
#             if random.random() > 0.5:
#                 print(1, end=" ", file=f)
#             else:
#                 print(0, end=" ", file=f)
#     print("", file=f)

print("----------graph generate finished!----------", file=sys.stderr)
print("----------graph generate finished!----------")
f.close()
