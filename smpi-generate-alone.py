#!/usr/bin/env python2.7
import sys
import os
import math

# Link parameters
link_latency = "10us"
link_bandwidth = 10
link_bandwidth_unit = "Gbps"

# XML generation functions
def issueHead():
        head = ("<?xml version='1.0'?>\n"
                "<!DOCTYPE platform SYSTEM \"http://simgrid.gforge.inria.fr/simgrid/simgrid.dtd\">\n"
                "<platform version=\"4.1\">\n\n")

        config_clause = ("<!--  WARNING:  This <config></config> clause below\n"
                       "makes it so that NO COMPUTATION TIME is simulated. This is because\n"
                       "in this module, for pedagogic purposes, we don't want to muddy the\n"
                       "(simulation) waters with computational times. As a results, this\n"
                       "XML platform file may not be suitable for running other\n"
                       "simulations, unless you remove the <config></config> clause.\n"
                       "-->\n"
                       "<config>\n"
                       "<prop id=\"smpi/simulate-computation\" value=\"0\"></prop>\n"
                       "<prop id=\"smpi/host-speed\" value=\""+str(real_compute_power)+"\"></prop>\n"
                       "</config>\n\n")

        AS_head = "<zone id=\"AS0\" routing=\"Full\">\n"

        return head + config_clause + AS_head

def issueHost(index):
  return "  <host id=\"host-"+str(index)+"."+hostname+"\" speed=\""+sim_compute_power+"\"/>\n"

def issueTail():
  return "</zone>\n</platform>\n"


hostname = "nimportequoi.fr"
real_compute_power = int(sys.argv[1]) #*1000000000
sim_compute_power = sys.argv[2]+"Gf"

filename = "./alone-platform.xml"
fh = open(filename, 'w')
fh.write(issueHead())

fh.write(issueHost(0))

fh.write(issueTail())
fh.close()

print >> sys.stderr, "a XML platform file created: "+filename

###############################################################
## Generate host file
filename = "./alone-hostfile.txt"
fh = open(filename, 'w')

fh.write("host-"+str(0)+"."+hostname+"\n")

fh.close()
print >> sys.stderr, "Hostfile created: "+filename


