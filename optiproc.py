import subprocess
import shlex
import matplotlib.pyplot as plt 
import numpy as np

def recupere_time(output):
    assert isinstance(output, str)

    lines = output.splitlines()
    #print(lines[-1])
    return float(lines[-1].split(" ")[-2])



def test(algo, nbprocs):
	#command = "mpirun -np "+str(nbprocs)+" ./graph-bisect  opti.txt bisect-"+algo+" 10 1.10 0"

	command = "smpirun -hostfile ring-"+str(nbprocs)+"-hostfile.txt -platform ring-"+str(nbprocs)+"-platform.xml ./graph  opti.txt bisect-"+algo+" 10 1.10 2"
	command = shlex.split(command)
	trial = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	output, stderr = trial.communicate()
	return recupere_time(str(output))


def create_ring(n):
	subprocess.Popen(["python", "smpi-generate-ring.py", str(n), "100", "100",  ".01Gbps", "1ms"])

nbprocs = []
for i in range(1,8):
	nbprocs.append(4*i)

for x in nbprocs:
	create_ring(x)

A = []
P = []
for i in range(7):
	print(i)
	A.append(test("a2a", nbprocs[i]))
	P.append(test("p2p", nbprocs[i]))

plt.figure(figsize=(10,5))
plt.title("Evolution with the number of processors for a ring")

plt.xlabel("Number of processors")
plt.ylabel("Time")
plt.plot(nbprocs, A, label="all-to-all", color='r')
plt.plot(nbprocs, P, label="peer-to-peer", color='g')
plt.legend()
plt.savefig("processors7.png")




