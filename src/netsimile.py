import time
import numpy as np
import scipy as sp
import matplotlib as mpl
import matplotlib.pyplot as plt
from pylab import *
import scipy.stats as st
import sys
from os import listdir
from os.path import isfile, join
from igraph import *
import scipy.spatial.distance
import matplotlib.pyplot as plt
import csv
from scipy.stats import entropy
from numpy.linalg import norm
import pandas as pd

#from utils import *

is_debug = False

# import local config: Set your local paths in dev_settings.py
DATA_URL=""
SAVE_URL=""
try:
    from dev_settings import *
except ImportError:
    pass


#==============================================================================
# 7 feature functions and their helper functions
# Page 2: Berlingerio, Michele, et al. "NetSimile: a scalable approach to
# size-independent network similarity." arXiv preprint arXiv:1209.2684 (2012).
#==============================================================================
def get_egonet(node, graph):
    """
    A nodes egonet is the induced subgraph formed by the node and its
    neighbors. Returns list of vertices that belong to the egonet
    """
    return graph.neighborhood(node)

def get_di(node, graph):
    """
    Number of neigbors
    """
    return graph.neighborhood_size(node)

def get_ci(node, graph):
    """
    Clustering coefficient of node, defined as the number of triangles
    connected to node over the number of connected triples centered on node
    """
    return graph.transitivity_local_undirected(node, mode=TRANSITIVITY_ZERO)

def get_dni(node, graph):
    """
    Average number of nodes two-hop away neighbors
    """
    return mean([get_di(n, graph) for n in get_egonet(node, graph) \
                                  if n != node])

def get_cni(node, graph):
    """
    Average clustering coefficient of neighbors of node
    """
    return mean([get_ci(n, graph) for n in get_egonet(node, graph) \
                                  if n != node])

def get_eegoi(node, graph):
    """
    Number of edges in nodes egonet;
    """
    edges_to=[]
    vertices = get_egonet(node, graph)

    for n in vertices:
        for i in get_egonet(n, graph):
            if i!=n:
                edges_to.append(i)

    #remove external nodes
    edges2x=[i for i in edges_to if i in vertices]
    assert (len(edges2x)%2==0),"Wrong calculation"
    return len(edges2x)/2


def get_eoegoi(node, graph):
    """
    Number of outgoing edges from node's egonet
    """
    edges_to=[]
    vertices=get_egonet(node, graph)

    for n in vertices:
        for i in get_egonet(n, graph):
            edges_to.append(i)

    return len([i for i in edges_to if i not in vertices])

def get_negoi(node, graph):
    """
    Number of neighbors of node's egonet
    """
    vertices = get_egonet(node, graph)
    all_neighbors = []
    for v in vertices:
        all_neighbors = all_neighbors + get_egonet(v,graph)
    all_neighbors = set(all_neighbors)
    all_neighbors =  [i for i in all_neighbors if i not in vertices]
    return(len(all_neighbors))

#==============================================================================
# KCore Algorithm components
#==============================================================================
def get_kshellProbs(g, name):
    print ("Extracting KCore features: %s" % name)
    #print(GraphBase.coreness(g))
    coreness = GraphBase.coreness(g)
    n = len(coreness)
    d = {n:coreness.count(n) for n in range(max(coreness)+1)} 
    kshellprobability = [d[key] / (n * 1.0) for key in sorted(d)]
    #print("kshellprob = " , kshellprobability)
    return (kshellprobability)


#==============================================================================
# NetSimile Algorithm components
# Features: Number of neigbors
#	    Clustering coefficient of node
#	    Average number of nodes two-hop away neighbors
# 	    Average clustering coefficient of neighbors of node
#	    Number of edges in nodes egonet
#	    Number of outgoing edges from node's egonet
#	    Number of neighbors of node's egonet
#==============================================================================
def get_features(g, name):
    print ("Extracting NetSimile features: %s" % name) 
    feature = [(get_di(i,g),
             get_ci(i,g),
             get_dni(i,g),
             get_cni(i,g),
             get_eegoi(i,g),
             get_eoegoi(i,g),
             get_negoi(i,g)) for i in g.vs]
    return feature

def get_features_all(graphs):
    """
    Returns all features of all graphs.
    Out Format: {g1:[(f1..f7),(f1..f7),(f1..f7)...#nodes], g2:...}
    """
    # Order all the graphs names based on the timestamp
    #ordered_names = sorted(graphs.keys(), key=lambda k:int(k.split('_',1)[0]))
    #return {g: get_features(graphs[g], g) for g in ordered_names}
    return {g: get_features(graphs[g], g) for g in graphs}

def get_moments(feat):
    """
    input: feature matrix of a single graph
    output: for each feature, return the 5 moments
    """
    #print("features: ", feat)
    feat_cols = zip(*feat)
    assert (len(feat_cols)==7),"Total columns != 7"

    # Calculate the 5 aggregates for each feature
    signature = []
    for f in feat_cols:
        print f
        signature = signature + [mean(f),
             median(f),
             std(f),
             st.skew(f),
             st.kurtosis(f)]
    return signature

def aggregator(features_all):
    #print("Aggregating features")
    return {g: get_moments(features_all[g]) for g in features_all}

def canberra_dist(sig1, sig2):
    """
    Returns the Canberra distance between graphs described by the signatures
    sig1 and sig2.
    """
    return abs(scipy.spatial.distance.canberra(sig1, sig2))

def KulbackLiebler_dist(sig1,sig2): #KulbackLiebler distance   sum(x *log(x/y))
    #maxLen = max(len(sig1), len(sig2))
    #sig1 = sig1 + [0]*(maxLen - len(sig1))
    #sig2 = sig2 + [0]*(maxLen - len(sig2))
    x = np.array(sig1)
    y = np.array(sig2)
    
    d = x*np.log2(x/y)
    d[np.isnan(d)] = 0
    d[np.isinf(d)] = 0
    #print 'KLD ', sum(d)
    return sum(d)

def JensenShannon_dist(sig1,sig2): #JensenShannon distance   sqrt(0.5 * KLD(x, (x+y)/2) + 0.5 * KLD(y, (x+y)/2))
    maxLen = max(len(sig1), len(sig2))
    sig1 = sig1 + [0]*(maxLen - len(sig1))
    sig2 = sig2 + [0]*(maxLen - len(sig2))
    x = np.array(sig1)
    y = np.array(sig2)
    m = (x+y) / 2
    d = sqrt(0.5 * KulbackLiebler_dist(x, m) + 0.5 * KulbackLiebler_dist(y, m))
    #print 'JSD', d
    return d
    

def saveDists(graph_names, dists, file_name):
    #assert (len(graph_names) - len(dists) == 1),\"len(graph_names) - len(dists) != 1"
    #print("graph names")
    #print(graph_names)
    data = zip(graph_names[0:],dists)
    #print("data: ", data)
    with open(file_name, 'w') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(graph_names)
        writer.writerows(data)
        #writer.writerows(dists)


def saveFeature(key, values, file_name):
    with open(file_name, 'a') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow([key] + values)
    return        

def saveFeatures(features, file_name):
    data = zip(features.keys(),features.values())
    with open(file_name, 'w') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerows(sorted(data))
                      
def find_NS_distance(sigs):
    sigkeys= sorted(sigs.keys())
    items = len(sigs) 
    out_slow = np.ones((items,items))
    for j in xrange(0, items):
        for k in xrange(j, items):
            out_slow[j, k] = canberra_dist(sigs[sigkeys[j]], sigs[sigkeys[k]]) 
            out_slow[k, j] = out_slow[j, k]
    return out_slow



def find_KC_distance(sigs):
    sigkeys= sorted(sigs.keys())
    items = len(sigs) 
    out_slow = np.ones((items,items))
    for j in xrange(0, items):
        for k in xrange(j+1, items):
            out_slow[j, k] = JensenShannon_dist(sigs[sigkeys[j]], sigs[sigkeys[k]]) 
            out_slow[k, j] = out_slow[j, k]
    return out_slow
     
def test_find_distance(sigs):
    sigkeys= sorted(sigs.keys())
    items = len(sigs) 
    out_slow = np.ones((items,items))
    for j in xrange(0,items):
        print("j = ", j , "key = ", sigkeys[j])
        for k in xrange(j, items):
            out_slow[j, k] = canberra_dist(sigs[sigkeys[j]], sigs[sigkeys[k]]) 
            out_slow[k, j] = out_slow[j, k]
            print("\tk = ", k , "key = ", sigkeys[k], "dist=", out_slow[j, k])
    print("out_slow = " , out_slow)
    return out_slow
                
def compareNS(sigs):
        for g in sigs:
            assert (len(sigs[g])==7*5),"Total features != 7*5"
            
        dists = find_NS_distance(sigs)
        #print("distance matrix");
        #print(dists)
        #print('ordered_graphs = ', ordered_graphs)
        saveDists(sorted(sigs.keys()), dists, sys.argv[1]+"_NSdists.txt")

def readDict(filename, sep):
    with open(filename, "r") as f:
        dict = {}
        for line in f:
            values = line.split(sep)
            dict[values[0]] = [float(x) for x in values[1:len(values)]]
        print dict
        return(dict)

def compareKCfromFile(filename):
    sigs = readDict(filename, ',')
    print("Read Signatures")
    #print(sigs)
    sigkeys= sigs.keys()
    l = len(sigkeys) 
    distanceMatrix = pd.DataFrame(np.ones((l,l)), index=sigkeys, columns=sigkeys)
    for i in range(l):
       g1name = sigkeys[i]
       for j in range(i,l):
            g2name = sigkeys[j]
            distanceMatrix[g1name][g2name] = jsd_dist(sigs[g1name], sigs[g2name]) 
            distanceMatrix[g2name][g1name] = distanceMatrix[g1name][g2name]
    distanceMatrix.to_csv(sys.argv[1].split('_')[0]+"_NSdists.txt")   
    return distanceMatrix        

def compareNSfromFile(filename):
    sigs = readDict(filename, ',')
    print("Read Signatures")
    #print(sigs)
    sigkeys= sigs.keys()
    l = len(sigkeys) 
    distanceMatrix = pd.DataFrame(np.ones((l,l)), index=sigkeys, columns=sigkeys)
    for i in range(l):
       g1name = sigkeys[i]
       for j in range(i,l):
            g2name = sigkeys[j]
            distanceMatrix[g1name][g2name] = canberra_dist(sigs[g1name], sigs[g2name]) 
            distanceMatrix[g2name][g1name] = distanceMatrix[g1name][g2name]
    distanceMatrix.to_csv(sys.argv[1].split('_')[0]+"_NSdists.txt")   
    return distanceMatrix
 
                 
def compareKC(sigs):           
        dists = find_KC_distance(sigs)
        #print("distance matrix");
        #print(dists)
        #print('ordered_graphs = ', ordered_graphs)
        saveDists(sorted(sigs.keys()), dists, sys.argv[1]+"_KCdists.txt")

def CompareNetSimileKCore(graph_files, dir_path):
    graphs = {f: Graph.Read_Ncol(join(dir_path, f), directed=False) for f in graph_files} 
    print("Read Graphs:",graphs.keys())
    
    NetSimileSignatures = {}
    KCoreSignatures = {}
    NetSimileVsKCoreTimings = {}
    
    ordered_graphs = sorted(graphs.keys())
    print("Ordered Graphs", ordered_graphs)
    
    t1 = 0
    t2 = 0
    for g in graphs:
       start_time = time.time()
       features = get_features(graphs[g], g)
       NetSimileSignature = get_moments(features)
       #print("NetSimile",g,NetSimileFeature)
       NetSimileSignatures[g] = NetSimileSignature
       saveFeature(g,NetSimileSignatures[g], sys.argv[1]+"_NetSimileSignatures.txt")
       t1 = time.time() - start_time
       
       start_time = time.time()
       KCoreSignature = get_kshellProbs(graphs[g], g)
       #print("KCoreAlgo",g, KCoreSignature)
       KCoreSignatures[g] = KCoreSignature
       t2 = time.time() - start_time
       NetSimileVsKCoreTimings[g] = [t1,t2]
    #saveFeatures(NetSimileSignatures, sys.argv[1]+"_NetSimileSignatures.txt")
    saveFeatures(KCoreSignatures, sys.argv[1]+"_KCoreSignatures.txt")
    #print (g,t1,'{0:f}'.format(t2), "seconds")
    saveFeatures(NetSimileVsKCoreTimings, sys.argv[1]+"_NetSimileVsKCoreTimings.txt")
    compareNS(NetSimileSignatures)
    compareKC(KCoreSignatures)
    return

def saveFeature(key, values, file_name):
    with open(file_name, 'a') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow([key] + values)
    return

def RuntimeNetSimileKCore(graph_files, dir_path):
    
    NetSimileVsKCoreTimings = {}   

    for g in graph_files:
       t0 = 0
       t1 = 0
       t2 = 0
       t3 = 0
       start_time = time.time()
       graph = Graph.Read_Ncol(join(dir_path, g))
       t0 = time.time() - start_time                     #time to read the graph
       
       print('Netsimile: Processing graph', g)
       start_time = time.time()
       features = get_features(graph, g)
       NetSimileSignature = get_moments(features)
       t1 = time.time() - start_time
       
       print('KCore: Processing graph', g)
       start_time = time.time()
       coreness = GraphBase.coreness(graph)
       n = len(coreness)
       d = {n:coreness.count(n) for n in range(max(coreness)+1)} 
       kshellprobability = [d[key] / (n * 1.0) for key in sorted(d)]
       t2 = time.time() - start_time
       
       start_time = time.time()
       maximumcore = max(coreness)
       EdgeMatrix = [[0 for x in range(maximumcore+1)] for x in range(maximumcore+1)] 
       for e in graph.es:
          i = coreness[e.source]
          j = coreness[e.target]
          EdgeMatrix[i][j] += 1
          if (i != j):
             EdgeMatrix[j][i] += 1
       numnodes = float(len(graph.vs))
       numedges = float(len(graph.es))
       Edgeprobabilities = []   
       if numedges == 0:
            numedges = 1
       for i in range(maximumcore+1):
          for j in range(i+1):
             Edgeprobabilities.append(EdgeMatrix[i][j]/numedges)

       t3 = time.time() - start_time
       
       NetSimileVsKCoreTimings[g] = [t0,t1,t2,t3]
       saveFeature(g, NetSimileVsKCoreTimings[g], sys.argv[1]+"_NetSimileVsKCoreTimings.txt")
    return



def getEdgeProbabilities(name, g, coreness, maximumcore):
    print ("Extracting EdgeProbabilities: %s" % name)
    EdgeMatrix = [[0 for x in range(maximumcore+1)] for x in range(maximumcore+1)] 
    for e in g.es:
       i = coreness[e.source]
       j = coreness[e.target]
       #print("Edge:", e.source,coreness[e.source],"-", e.target,coreness[e.target])
       EdgeMatrix[i][j] += 1
       if (i != j):
         EdgeMatrix[j][i] += 1
    
    numnodes = float(len(g.vs))
    numedges = float(len(g.es))
    #print("nodes:", numnodes, " edges: ", numedges)
    #for i in range(maximumcore+1):
    #   for j in range(maximumcore+1):
    #      EdgeMatrix[i][j] = EdgeMatrix[i][j]/numedges
 
    Edgeprobabilities = []   
    if numedges == 0:
      numedges = 1
    for i in range(maximumcore+1):
       for j in range(i+1):
          Edgeprobabilities.append(EdgeMatrix[i][j]/numedges)
          
    #print ("EdgeMatrix:", np.tril(EdgeMatrix))
    #return np.tril(EdgeMatrix).tolist()
    #print ("Edgeprobabilities", Edgeprobabilities)
    return Edgeprobabilities
#==============================================================================
# NetSimile algorithm
#==============================================================================
def NetSimile(graph_files, dir_path, use_old_dists=False):
    
    NetSimileSignatures = {}
    NetSimileTimings = {}

    for g in graph_files:
       graph = Graph.Read_Ncol(join(dir_path, g), directed=False)
       start_time = time.time()
       features = get_features(graph, g)
       NetSimileSignature = get_moments(features)
       NetSimileSignatures[g] = NetSimileSignature      
       NetSimileTimings[g] = time.time() - start_time
       saveFeature(g, NetSimileSignature, sys.argv[1]+"_NetSimileSignatures.txt")
    saveFeatures(NetSimileTimings, sys.argv[1]+"_NetSimileTimings.txt")
    dists = find_NS_distance(NetSimileSignatures)
    saveDists(sorted(NetSimileSignatures.keys()), dists, sys.argv[1]+"_NSdists.txt")
    return

#==============================================================================
# Main
# Command line parameter: name-of-dataset
# Example Usage:$ python netsimile.py "reality_mining_voices"
#==============================================================================
if __name__=="__main__":
    #dir_path = join(DATA_URL, sys.argv[1])
    #print(dir_path)
    #onlyfiles = [f for f in listdir(dir_path) if isfile(join(dir_path,f)) ]
    #print("Files in Directory: ",onlyfiles)
    #NetSimile(onlyfiles, dir_path, use_old_dists=is_debug)
    #CompareNetSimileKCore(onlyfiles, dir_path)
    compareNSfromFile(sys.argv[1])
    #RuntimeNetSimileKCore(onlyfiles, dir_path)
    #filename = sys.argv[1] + "NS-Data_NetSimileSignatures.txt"
    #print filename
    #compareNSfromFile(filename)
    

