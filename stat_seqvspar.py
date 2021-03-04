import subprocess
import shlex
import matplotlib.pyplot as plt 
import numpy as np

def recupere_time(output):
    assert isinstance(output, str)

    lines = output.splitlines()
    #print(lines[-1])
    return float(lines[-1].split(" ")[-2])


def create_graphe(n,m):
	name = "stats/graph_full"+str(n)+"-"+str(m)+".txt"
	subprocess.Popen(["./create-graph.py",str(n), str(m),name])

def create_all(N,E):
	for i in range(len(N)):
		create_graphe(N[i], E[i])
	print("graphs succesfully created")

def test(n,m, algo, nbprocs):
	command = "mpirun -np "+str(nbprocs)+" ./graph-bisect stats/graph_full"+str(n)+"-"+str(m)+".txt bisect-"+algo+" 10 1.10"
	command = shlex.split(command)
	trial = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	output, stderr = trial.communicate()
	return recupere_time(str(output))


"""
Premiere partie : creation des graphs
avec un deux fois plus d arretes que de sommets

"""

N = [100*i for i in range(1,40, 2)]
#E = [2*N[i] for i in range(len(N))]
#E = [int(2*np.log(N[i])) for i in range(len(N))]
E = [int(N[i]*N[i]/4) for i in range(len(N))]

create_all(N,E)


# Deuxieme partie : tourner sur les trois algos

t_seq=[]
t_all= []
t_peer=[]

for i in range(len(N)):
	t_seq.append(test(N[i],E[i], "seq", 1))
	t_all.append(test(N[i], E[i], "a2a", 4))
	t_peer.append(test(N[i], E[i], "p2p", 4))


plt.figure(figsize=(10,5))
#plt.title("Comparison for graph with a linear number of edges using mpirun")
#plt.title("Comparison for graph with a logarithmic number of edges using mpirun")
plt.title("Comparison for graph with a quadratic number of edges using mpirun")

plt.xlabel("Number of vertices")
plt.ylabel("Time")
plt.plot(N, t_seq, label="seq", color='b')
plt.plot(N, t_all, label="all-to-all", color='r')
plt.plot(N, t_peer, label="peer-to-peer", color='g')
plt.legend()
plt.savefig("full_graphs.png")




