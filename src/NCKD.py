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
import matplotlib.cm as cm
import csv
from scipy.stats import entropy
from numpy.linalg import norm
#from cv2 import *
import pandas as pd

from utils import *

# import local config: Set your local paths in dev_settings.py
DATA_URL=""
SAVE_URL=""
try:
    from dev_settings import *
except ImportError:
    pass
    
#==============================================================================
# KCore Algorithm components
#==============================================================================

    
def get_kshellProbs(g, name):
    print ("Extracting KCore features: %s" % name)
    #print(GraphBase.coreness(g))
    coreness = GraphBase.coreness(g)
    n = len(coreness)
    d = {n:coreness.count(n) for n in range(max(coreness)+1)} 
    #print("node distn = " , d)
    kshellprobability = [d[key] / (n * 1.0) for key in sorted(d)]
    #print("kshellprob = " , kshellprobability)
    return (kshellprobability)
   
def getEdgeMatrix(g, name):
    print ("Extracting EdgeMatrix: %s" % name)
    #print("Coreness: ",GraphBase.coreness(g))
    coreness = GraphBase.coreness(g)
    highestcore = max(coreness)
    g.vs["coreness"] = coreness
    EdgeMatrix = [[0 for x in range(highestcore+1)] for x in range(highestcore+1)] 
    for e in g.es:
       i = coreness[e.source]
       j = coreness[e.target]
       #print("Edge:", e.source,coreness[e.source],"-", e.target,coreness[e.target])
       EdgeMatrix[i][j] += 1
       if (i != j):
         EdgeMatrix[j][i] += 1
    print EdgeMatrix
    return EdgeMatrix

    
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


def NodeDistribution(graphs):
    KCoreSignatures = {}
    KCoreTimings = {}   
    t1 = 0
    for g in graphs:
       start_time = time.time()
       KCoreSignature = get_kshellProbs(graphs[g], g)
       KCoreSignatures[g] = KCoreSignature
       t1 = time.time() - start_time
       KCoreTimings[g] = [t1]
       saveFeature(g,KCoreSignatures[g], sys.argv[1]+"_NodeDistn.txt")
    #print (g,t1,'{0:f}'.format(t2), "seconds")
    saveDict(KCoreTimings, sys.argv[1]+"_NodeDistnTimings.txt")
    return
    
def NodeDistributionGiant(graphs):  
    KCoreSignatures = {}
    KCoreTimings = {}   
    t1 = 0
    for g in graphs:
       G = graphs[g].components().giant()
       start_time = time.time()
       KCoreSignature = get_kshellProbs(G, g)
       #print("KCoreAlgo",g, KCoreSignature)
       KCoreSignatures[g] = KCoreSignature
       t1 = time.time() - start_time
       KCoreTimings[g] = [t1]
    saveDict(KCoreSignatures, sys.argv[1]+"_Giant_NodeDistn.txt")
    #print (g,t1,'{0:f}'.format(t2), "seconds")
    saveDict(KCoreTimings, sys.argv[1]+"_Giant_NodeDistnTimings.txt")
    return

        
def ComputeEdgeMatrix(graphs):
    Coreness = {}
    EdgeMatrices = {}
    TimeToCompute = {}
   
    maximumcore = 0
    for g in graphs:
       Coreness[g] = GraphBase.coreness(graphs[g])
       highestcore = max(Coreness[g])
       if (maximumcore < highestcore):
           maximumcore = highestcore    
    print("maximumcore", maximumcore)
    
    for g in graphs:
       #plotKCore(graphs[g])
       start_time = time.time()
       EdgeMatrices[g] = getEdgeProbabilities(g, graphs[g], Coreness[g], maximumcore)
       TimeToCompute[g] = time.time() - start_time
       #get_kshellProbs(graphs[g],g)
       saveFeature(g,EdgeMatrices[g], sys.argv[1]+"_EdgeProbabilities.txt")
    saveDict(TimeToCompute, sys.argv[1]+"_EdgeProbabilitiesTiming.txt")
    #print("EdgeMatrices: " , EdgeMatrices)
    #print("Timing: ", TimeToCompute)
    
    return


def ComputeEdgeMatrixGiant(graphs):
    Coreness = {}
    EdgeMatrices = {}
    TimeToCompute = {}
     
    maximumcore = 0
    for g in graphs:
       G = graphs[g].components().giant()
       Coreness[g] = getCoreness(G, g)
       highestcore = max(Coreness[g])
       if (maximumcore < highestcore):
           maximumcore = highestcore    
    print("maximumcore", maximumcore)
    
    for g in graphs:
       start_time = time.time()
       G = graphs[g].components().giant()
       EdgeMatrices[g] = getEdgeProbabilities(g, G , Coreness[g], maximumcore)
       TimeToCompute[g] = time.time() - start_time
       
    saveDict(EdgeMatrices, sys.argv[1]+"_Giant_EdgeProbabilities.txt")
    saveDict(TimeToCompute, sys.argv[1]+"_Giant_EdgeProbabilitiesTiming.txt")
    #print("EdgeMatrices: " , EdgeMatrices)
    #print("Timing: ", TimeToCompute)

    return

def getComponentProbabilities(name, g, coreness):
    print ("Extracting Component Probabilities: %s" % name)
    ComponentProbabilities = []
    highestcore = max(coreness)
    kshell = [[] for k in range(highestcore+1)] 
    for v in g.vs:
      i = coreness[v.index]
      kshell[i].append(v)
    
    NumberOfComponents = [0 for k in range(highestcore+1)]
    ComponentSizes = [{} for k in range(highestcore+1)]
    ShellFeatures = [{} for k in range(highestcore+1)]
    
    for k in range(highestcore+1):
      subgraph = g.subgraph(kshell[k])
      componentslist = subgraph.components()
      NumberOfComponents[k] = len(componentslist)
      ShellFeatures[k] = [subgraph.average_path_length(), subgraph.diameter(), subgraph.radius(), subgraph.transitivity_undirected()]
      cs = {}
      #cs[0] =  len(componentslist)
      for component in componentslist:
         size = len(component)
         if not size in cs:
           cs[size] = 1
         else:
           cs[size] += 1
           
      ComponentSizes[k] = cs
      
      
      #print("k = ", k)     
      #print("subgraph vertices")
      #for v in subgraph.vs:
      #   print(v['name']);
      #print("subgraph edges")
      #for e in subgraph.es:
      #  print (subgraph.vs[e.source]['name'], "-- ", subgraph.vs[e.target]['name'])
      #componentslist = subgraph.components()
      #print("num components = ", len(componentslist))
      #for component in componentslist:
      #   print component
    #print ComponentSizes
    #print ShellFeatures
    saveFeature(name,NumberOfComponents, sys.argv[1]+"_NumberOfComponents.txt")
    saveFeature(name,ComponentSizes, sys.argv[1]+"_ComponentSizes.txt")
    saveFeature(name,ShellFeatures, sys.argv[1]+"_ShellFeatures.txt") 
    return NumberOfComponents,ComponentSizes, ShellFeatures
    

def ComputeComponentprobabilities(graphs):
    ComponentSizes = {}
    NumberOfComponents = {}
    ShellFeatures = {}
    
    for g in graphs:
       NumberOfComponents[g], ComponentSizes[g], ShellFeatures[g] = getComponentProbabilities(g, graphs[g], GraphBase.coreness(graphs[g]))
 
    #saveDict(NumberOfComponents, sys.argv[1]+"_NumberOfComponents.txt")
    #saveDict(ComponentSizes, sys.argv[1]+"_ComponentSizes.txt")
    #saveDict(ShellFeatures, sys.argv[1]+"_ShellFeatures.txt")   
    return



def NodeEdgeDistribution(graphs):
    #graphs = {}
    Coreness = {}
    KCoreSignatures = {}
    KCoreTimings = {}
    CountTimings = {}
    NodeDistnTimings = {}   
    EdgeMatrices = {}
    EdgeDistnTimings = {}
    AllTimings = {}
  
    nodedistnfilename = sys.argv[1]+"_NodeDistnProbabilities.txt"
    edgedistnfilename = sys.argv[1]+"_EdgeProbabilities.txt"
    os.remove(nodedistnfilename) if os.path.exists(nodedistnfilename) else None
    os.remove(edgedistnfilename) if os.path.exists(edgedistnfilename) else None
    
    maximumcore = 0
    for g in graphs:
       print('Processing graph:', g)
       #graph = Graph.Read_Ncol(join(dir_path, g), directed=False)  
       print ("Extracting Node Distn Probabilities: %s" % g)
       graph = graphs[g].simplify()
       start_time = time.time()
       coreness = GraphBase.coreness(graph)
       KCoreTimings[g] = time.time() - start_time
       Coreness[g] = coreness
       
       l = len(coreness)
       
       start_time = time.time()
       d = {n:coreness.count(n) for n in range(max(coreness)+1)}
       CountTimings[g] = time.time() - start_time
       print('max(coreness): ', max(coreness))
       start_time = time.time()
       KCoreSignature = [d[key] / (l * 1.0) for key in sorted(d)]
       NodeDistnTimings[g] = time.time() - start_time
       
       KCoreSignatures[g] = KCoreSignature
       
       saveFeature(g,KCoreSignature, nodedistnfilename)
       highestcore = max(Coreness[g])
       if (maximumcore < highestcore):
           maximumcore = highestcore    
       
    saveDict(NodeDistnTimings, sys.argv[1]+"_NodeDistnTimings.txt") 
          
 
    for g in graphs:
       start_time = time.time()
       EdgeMatrix = getEdgeProbabilities(g, graphs[g], Coreness[g], maximumcore)
       EdgeMatrices[g] = EdgeMatrix
       EdgeDistnTimings[g] = time.time() - start_time
       saveFeature(g,EdgeMatrix, edgedistnfilename)
    saveDict(EdgeDistnTimings, sys.argv[1]+"_EdgeProbabilitiesTiming.txt")    
    
    for g in graphs:
       t0= KCoreTimings[g]
       t1 = CountTimings[g]
       t2 = NodeDistnTimings[g]
       t3 = EdgeDistnTimings[g]
       
       t4 =  t0+t1+t2+t3
       AllTimings[g]=[t0,t1 ,t2, t3, t4]   
       #AllTimings[g]=[t3]
    
    distance(nodedistnfilename)
    distance(edgedistnfilename)
    saveDict(AllTimings, sys.argv[1]+"_AllNCKDTimings.txt")  
    return

def NewNodeEdgeDistibutions(graphs):   # node + intra-edge
    nodedistnfilename = sys.argv[1]+"_NewNodeDistnProbabilities.txt"
    intraedgedistnfilename = sys.argv[1]+"_NewIntraEdgeDistnProbabilities.txt"
    os.remove(nodedistnfilename) if os.path.exists(nodedistnfilename) else None
    os.remove(intraedgedistnfilename) if os.path.exists(intraedgedistnfilename) else None
    
    computationtime = {}
    for g in graphs:
       print('Processing graph:', g)
       #graph = Graph.Read_Ncol(join(dir_path, g), directed=False)  
       print ("Extracting Node _ IntraEdge Distn Probabilities: %s" % g)
       graph = graphs[g] #.simplify()
       start_time = time.time()
       coreness = GraphBase.coreness(graph)
       highestcore = max(coreness)
       kshell = [[] for k in range(highestcore+1)] 
       for v in graph.vs:
          i = coreness[v.index]
          kshell[i].append(v) 
       
       nodedistn = []
       intraedgedistn = []
       
       for k in range(highestcore+1):
          subgraph = graph.subgraph(kshell[k])
          subgraph.simplify()
          nodedistn.append(float(subgraph.vcount()))
          intraedgedistn.append(float(subgraph.ecount()))
       
       sumn = sum(nodedistn)
       normalizednodedistn = []   
       if not(sumn ==0):
           normalizednodedistn = [x/sumn for x in nodedistn]

       sumn = sum(intraedgedistn)
       normalizedintraedgedistn = []   
       if not(sumn ==0):
           normalizedintraedgedistn = [x/sumn for x in intraedgedistn]
           
       computationtime[g] = time.time() - start_time
      
       saveFeature(g,normalizednodedistn, nodedistnfilename)
       saveFeature(g,normalizedintraedgedistn, intraedgedistnfilename)

    distance(nodedistnfilename)
    distance(intraedgedistnfilename)
    saveDict(computationtime, sys.argv[1]+"_NewNCKDTimings.txt")   
    return
    
def distance(filename):
    outputfile = os.path.splitext(filename)[0] +"_JSDistanceMatrix.txt" 
    print "outputfle ", outputfile
    os.remove(outputfile) if os.path.exists(outputfile) else None
    
    sigs = readDict(filename, ',')
    print "Read Distributions", filename
    sigkeys= sigs.keys()
    l = len(sigkeys) 
    JSDistanceMatrix = pd.DataFrame(np.ones((l,l)), index=sigkeys, columns=sigkeys)
    for i in range(l):
       g1name = sigkeys[i]
       for j in range(i,l):
            g2name = sigkeys[j]
            JSDistanceMatrix[g1name][g2name] = JensenShannon_dist(sigs[g1name], sigs[g2name]) 
            JSDistanceMatrix[g2name][g1name] = JSDistanceMatrix[g1name][g2name] 
    JSDistanceMatrix.to_csv(outputfile)   
    return JSDistanceMatrix    
       

def Runtime(dir_path,graphfiles):
    #graphs = {}
    Coreness = {}
    KCoreSignatures = {}
    KCoreTimings = {}
    CountTimings = {}
    NodeDistnTimings = {}   
    EdgeMatrices = {}
    EdgeDistnTimings = {}
    AllTimings = {}

    for g in graphfiles:
       print('Processing graph:', g)
       graph = Graph.Read_Ncol(join(dir_path, g), directed=False)  
       print ("Extracting Node Distn Probabilities: %s" % g)
       start_time = time.time()
       coreness = GraphBase.coreness(graph)
       Coreness[g] = coreness
       l = len(coreness)
       d = {n:coreness.count(n) for n in range(max(coreness)+1)}
       KCoreSignature = [d[key] / (l * 1.0) for key in sorted(d)]
       NodeDistnTimings[g] = time.time() - start_time
       KCoreSignatures[g] = KCoreSignature
       print('max(coreness): ', max(coreness))
       saveFeature(g,KCoreSignature, sys.argv[1]+"_NodeDistnProbabilities.txt")
       start_time = time.time()
       EdgeMatrix = getEdgeProbabilities(g, graph, Coreness[g], max(coreness))
       EdgeMatrices[g] = EdgeMatrix
       EdgeDistnTimings[g] = time.time() - start_time
       saveFeature(g,EdgeMatrix, sys.argv[1]+"_EdgeProbabilities.txt")
       AllTimings[g]= NodeDistnTimings[g] + EdgeDistnTimings[g]
       
    saveDict(NodeDistnTimings, sys.argv[1]+"_NodeDistnTimings.txt") 
    saveDict(EdgeDistnTimings, sys.argv[1]+"_EdgeProbabilitiesTiming.txt")  
    saveDict(AllTimings, sys.argv[1]+"_AllNCKDTimings.txt")        
    return



def readDict(filename, sep):
    with open(filename, "r") as f:
        dict = {}
        for line in f:
            values = line.split(sep)
            dict[values[0]] = [float(x) for x in values[1:len(values)]]
        return(dict)

def saveDists(graph_names, dists, file_name):
    print('data originally')
    print(dists)    
    
    print('data as list')
    #print(list(dists))
    
    print(data)
    data = zip(graph_names[0:],dists)
    with open(file_name, 'w') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(graph_names)
        writer.writerows(data)
    return
                
def emd(x,y):
  #print("x = ", x ," y = " , y) 
  a = np.zeros((len(x),2))
  b = np.zeros((len(y),2))
  for i in range(0,len(a)):
    a[i][0] = x[i]
    a[i][1] = i+1.0

  for i in range(0,len(b)):
    b[i][0] = y[i]
    b[i][1] = i+1.0
  #print("a = ", a ," b = " , b) 
  # Convert from numpy array to CV_32FC1 Mat
  a64 = cv.fromarray(a)
  a32 = cv.CreateMat(a64.rows, a64.cols, cv.CV_32FC1)
  cv.Convert(a64, a32)

  b64 = cv.fromarray(b)
  b32 = cv.CreateMat(b64.rows, b64.cols, cv.CV_32FC1)
  cv.Convert(b64, b32)

  # Calculate Earth Mover's
  #print cv.CalcEMD2(a32,b32,cv.CV_DIST_L2)
  return cv.CalcEMD2(a32,b32,cv.CV_DIST_L2)
  
  
def find_EMD_distance(filename):
    sigs = readDict(filename, ',')
    sigkeys= sorted(sigs.keys())
    items = len(sigs) 
    dists = np.ones((items,items))
    for j in xrange(0, items):
        for k in xrange(j, items):
            a = np.array(sigs[sigkeys[j]])
            b = np.array(sigs[sigkeys[k]])
            dists[j, k] = emd(a, b) 
            dists[k, j] = dists[j, k]
    print(dists)     
    saveDists(sorted(sigs.keys()), dists, sys.argv[1].split('_')[0]+"_NodeEMDdists.txt")
    return dists
    
def new_find_EMD_distance(filename):
    sigs = readDict(filename, ',')
    print("Read Signatures")
    print(sigs)
    sigkeys= sigs.keys()
    l = len(sigkeys) 
    distanceMatrix = pd.DataFrame(np.ones((l,l)), index=sigkeys, columns=sigkeys)
    for i in range(l):
       g1name = sigkeys[i]
       for j in range(i,l):
            g2name = sigkeys[j]
            print('g1:', g1name)
            print('g2:', g2name)
            a = np.array(sigs[g1name])
            b = np.array(sigs[g2name])
            print('a:' , a)
            print('b:', b)
            print('emd', emd(a,b))
            distanceMatrix[g1name][g2name] = emd(a, b) 
            distanceMatrix[g2name][g1name] = distanceMatrix[g1name][g2name]
    distanceMatrix.to_csv(sys.argv[1].split('_')[0]+"_NodeEMDdists.txt")   
    return distanceMatrix
     

def Test():
    g1 = Graph.Full(20)
    g2 = Graph.Full(20)
    g2.delete_edges(g2.es)
    
    coreness1 = GraphBase.coreness(g1)
    coreness2 = GraphBase.coreness(g2)
    
    n1 = len(coreness1)
    d1 = {n:coreness1.count(n) for n in range(max(coreness1)+1)} 
    KCoreSignature1 = [d1[key] / (n1 * 1.0) for key in sorted(d1)]    
    n2 = len(coreness2)
    d2 = {n:coreness2.count(n) for n in range(max(coreness2)+1)} 
    KCoreSignature2 = [d2[key] / (n2 * 1.0) for key in sorted(d2)]   
     
    print "KCoreSignature1:", KCoreSignature1
    print "KCoreSignature2: ", KCoreSignature2
    
    maximumcore = max(max(coreness1),max(coreness2))
    print "maximumcore", maximumcore
    EdgeMatrix1 = getEdgeProbabilities("full", g1, coreness1, maximumcore)
    EdgeMatrix2 = getEdgeProbabilities("empty", g2, coreness2, maximumcore)
    
    print "EdgeMatrix1", EdgeMatrix1
    print "EdgeMatrix2", EdgeMatrix2
    return   
    
def testSingle():
       g = sys.argv[1]
       print "Processing graph:", g
       graph = Graph.Read_Ncol(g, directed=False)
       print "Extracting Node Distn Probabilities: %s" % g
       coreness = GraphBase.coreness(graph)
       l = len(coreness)
       d = {n:coreness.count(n) for n in range(max(coreness)+1)}
       NodeProbabilities = [d[key] / (l * 1.0) for key in sorted(d)]

       print 'NCKD NodeProbabilities'
       print  NodeProbabilities
       
       EdgeProbabilities = getEdgeProbabilities(g, graph, coreness, max(coreness))
       print 'NCKD EdgeProbabilities'
       print  EdgeProbabilities
       
       NumberOfComponents,ComponentSizes, ShellFeatures = getComponentProbabilities(g, graph, coreness)
       print 'NCKD ComponentProbabilities'
       print  'NumComponents', NumberOfComponents
       print 'ComponentSizes', ComponentSizes
       print 'ShellFeatures', ShellFeatures
       return 

def runDir():
       dir_path = join(DATA_URL, sys.argv[1])
       print(dir_path)
       graph_files = [f for f in listdir(dir_path) if \
                            isfile(join(dir_path,f)) ]
       graphs = {f: Graph.Read_Ncol(join(dir_path, f), directed=False) for f in graph_files} 
       print("Read Graphs:",graphs.keys())
       #NewNodeEdgeDistibutions(graphs)
       NodeEdgeDistribution(graphs)
       #Runtime(dir_path, graph_files)
       #NodeDistribution(graphs)
       #ComputeEdgeMatrix(graphs)
       #ComputeComponentprobabilities(graphs)
       return 
    
#==============================================================================
# Main
# Command line parameter: name-of-dataset
# Example Usage:$ python NCKD.py "./TestData/data"
#==============================================================================
if __name__=="__main__":
    #Test()
    #testSingle()
    runDir()
    #new_find_EMD_distance(sys.argv[1])
    #nodedistnfilename = sys.argv[1] + "NS-Data_NodeDistnProbabilities.txt"
    #print nodedistnfilename
    #distance(sys.argv[1])
