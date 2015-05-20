
import sys
from pathway import *
from geneexp import *
from graphs import *
from os import listdir
import networkx as nx
import random
import matplotlib.pyplot as plt
from fdrcorrection import fdrcorrection0
import numpy as np

expf = sys.argv[1]; # expression file

genenameidx = 2; # default value = 0 (set it to 0 if user doesn't input any val)
foldchangeidx = 9; # default value = 1 # IF LOG FOLD CHANGE, ELIMINATE INF AND -INF
state1idx = 7; # NI
state2idx = 8; # ND # ALLOW USER TO INPUT DATA AS EXPRESSION VALUE OR RATIO
nrand = 1000; # number of randomizations for calculating p-values

nodelabels = ['p','pp','ppp','n','nn','nnn'];
edgelabels = ['activate','inhibit','interact'];
nekeys = [];
connectivitydict = {};
idx = 0;
for e in edgelabels:
    for n1 in nodelabels:
        for n2 in nodelabels:
            key = n1+"-"+n2+"-"+e;
            if (n1 == n2): # unique index
                connectivitydict[key] = idx;
                idx += 1;
            else: # if labels are diff, then get same id for symmetric labels, i.e. p-n-activate gets same label as n-p-activate
                symmetrickey = n2+"-"+n1+"-"+e;
                if symmetrickey in connectivitydict:
                    connectivitydict[key] = connectivitydict[symmetrickey]; # get index of symmetric key
                else:
                    connectivitydict[key] = idx;
                    idx += 1;

exp_obj = geneexp(expf, genenameidx, foldchangeidx, state1idx, state2idx); 
genelabels = exp_obj.assignNodeLabels(); # assign labels to gene based on gene expression
print exp_obj.labeldist;
labelprob = {k : float(exp_obj.labeldist[k])/sum(exp_obj.labeldist.values()) for k in exp_obj.labeldist.keys()}
print labelprob;
# print genelabels

def expLabelsPathgenes( pgenedict ):
    pgenelabels = {}; # this dict will be used to relabel pathway graph with gene-exp labels instead of entry-ids from kegg xml files
#    print pgenedict

    for id in pgenedict: # id is key, that is used as node label in the pathway graph. value is gene name/aliases
        pgenenames = pgenedict[id].split(', ');
        if ('...' in pgenenames[-1]):
            del pgenenames[-1]; # sometimes gene name list ends with '...', and it seems the last gene alias is not complete when it ends with ..., so remove it.

        # find gene names in genelabels and get the label(s)
        pglabel = [];
        for egenes in genelabels:
            eglist = egenes.split(','); # genes in the exp file can also be in a comma-separated list
            for eg in eglist:
                for pg in pgenenames:
                    if (eg == pg): # only if pathway gene (pg) is an exact match to the exp-list gene (eg)
                        pglabel.append( genelabels[egenes] );
                        continue;
#        print id, pgenenames, pglabel
        # if pathway genes have no label (i.e. no expression data) or multiple labels, remove those genes (nodes) from the pathway before generating subgraphs        
        if (len(pglabel) == 1): # each node/pgene should have exactly one label
            pgenelabels[id] = pglabel[0]; # assign gene exp label to entry id (node id)

    return pgenelabels;
	
def randomizeLabels( pg ): # read a dict with key: node is, and value: exp label, and return new dict with randomized labels
    # randomize the mapping between gene ids and expression labels: pick labels from the global label distribution
    ngenes = len(pg);
    randlabels = ['p'] * int(round(labelprob['p'] * ngenes)) + ['pp'] * int(round(labelprob['pp'] * ngenes)) + ['ppp'] * int(round(labelprob['ppp'] * ngenes)) + ['n'] * int(round(labelprob['n'] * ngenes)) + ['nn'] * int(round(labelprob['nn'] * ngenes)) + ['nnn'] * int(round(labelprob['nnn'] * ngenes));
    lendiff = len(pg) - len(randlabels);
    if (lendiff > 0): # fewer labels in randlabels. 'n' or 'p' should be most common: add that many n/p labels
        randlabels = randlabels + ['n'] * lendiff;
    elif (lendiff < 0): # more labels in randlabels, remove from end
        del randlabels[lendiff:]; # removes last absdiff elements from the list
    random.shuffle( randlabels );
    return dict( zip( pg.keys(), randlabels ) ); # randomize keys in pg and store that as a new dict

def retrievePathway( listsubg, genedict ):
    for g in listsubg:
        # print g.edges(data=True);
        ginfo = g.edges(data=True); # list of edges in this subgraph 
        path = [];
        for c in ginfo: # list where each item is (n1, n2, {'type': 'edgetype'})
            path.append( "["+genedict[c[0]]+"] <--"+c[2]['type']+"--> ["+genedict[c[1]]+"]" );
        return ";".join( path );
        #print genedict[ ginfo[0] ], genedict[ ginfo[1] ], ginfo[2]['type']
 
numsubgraphs = 0; # number of subgraphs in all pathways
allkgmlfiles = [f for f in listdir("./KEGGpathways/")]

pvals = [];
sigsubgraphs = {}; # key is pathway name, val is dict: key is p-value, and value is list of subgraphs with that p-value

for f in allkgmlfiles:
    if (".xml" not in f):
        continue  # only read .xml files
    pobj = pathway("./KEGGpathways/"+f);
    print "\n===\n",f, pobj.pathname;
#    print "genedict: ", pobj.genedict;
    pgenes = expLabelsPathgenes( pobj.genedict ); # pgenes has node id as key and exp label as value
#    print "mapping: ", pgenes;
    pobj.removeNodes( pgenes ); # remove all nodes that do not have a mapping to a gene exp label
#    print pobj.pathwaygraph.edges(data=True);
    #pos = nx.spring_layout(pobj.pathwaygraph);
    #nx.draw(pobj.pathwaygraph, pos);
    #plt.savefig("graph.pdf"); 
    #plt.close();
    psubgraphs = pobj.getSubgraphs( 4 ); # generate all connected subgraphs of size 3
    if ( len(psubgraphs) == 0 ):
        print "No connected subgraphs in this pathway.\n";
        continue;
    numsubgraphs += len(psubgraphs);
    gdb = graphs(); # create new graph instance for this pathway
    pgidmap = gdb.populateGDB( psubgraphs, pgenes, connectivitydict ); # pgidmap:: key: subgraph gid, value: list of subgraphs with that gid
    gdb.getStats();

    # initialize the prandcounts (how many times each gid is observed when node labels are randomized)
    prandcounts = {k : 0 for k in pgidmap.keys()};        
    sigsubgraphs[pobj.pathname] = {};

    print "Randomizing gene expression labels to compute p-values ...\n";
    for i in range(0,nrand):
        # create a new instance of gdb for this randomization
        gdbrand = graphs();
        randpg = randomizeLabels( pgenes ); # randomize gene exp labels for genes in this pathway
        randgidmap = gdbrand.populateGDB( psubgraphs, randpg, connectivitydict ); # provide subgraphs and randomized gene-expression mapping

        for gid in pgidmap:
            if gid in gdbrand.gdb:
                # If a graph occurs in randomization, count it only once instead of the number of times it occurs
                prandcounts[gid] += 1; # keeping track of graph occuring in all randomizations
      
#    print "From ", nrand, " randomizations of gene expression labels:";
    for gid in prandcounts:
#        print gid, len(pgidmap[gid]), prandcounts[gid];
        pvalue = float(prandcounts[gid])/nrand;
        pvals.append( pvalue );
        
        # if p-value < 0.05, retrieve the corresponding subgraph
        if (pvalue < 0.05):
            sigpath = retrievePathway( pgidmap[gid], pobj.genedict );
            # print sigpath;
            if pvalue in sigsubgraphs[pobj.pathname]:
                sigsubgraphs[pobj.pathname][pvalue].append( sigpath );
            else:
                sigsubgraphs[pobj.pathname][pvalue] = [sigpath];

binaryarr, pvals_corr = fdrcorrection0(pvals, 0.1, method='indep'); # 10% FDR rate (GSEA uses 25%)
#print pvals;
#print binaryarr;
#print pvals_corr;

plt.hist(pvals);
plt.savefig("pval_hist.pdf");
plt.close();

plt.hist(pvals_corr);
plt.savefig("pval_corr_hist.pdf");
plt.close();

numsigpath = 0;
for pname in sigsubgraphs:
    sigkeys = sigsubgraphs[pname].keys();
    print "\n", pname;
    for pv in sigsubgraphs[pname].keys():
        if (pv < 0.009):
            numsigpath += 1;
            print "p-value: ", pv;
            for sigp in sigsubgraphs[pname][pv]:
                print sigp, "\n";
#    print sigsubgraphs[pname];
print "Number of significant pathways (original p-value < 0.009):", numsigpath;
print "Total number of subgraphs analyzed: ", numsubgraphs;
