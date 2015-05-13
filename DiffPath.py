
import sys
from pathway import *
from geneexp import *
from graphs import *
from os import listdir
import networkx as nx
import pickle

#pathwayf = sys.argv[1];
expf = sys.argv[1]; # expression file

#path_obj = pathway(pathwayf); # read kgml pathway file
#obj.getSubPathways(); # get connected components in the pathway file
#sub3 = path_obj.getSubgraphs( 3 );
#print sub3;
#pickle.dump(sub3, open('sub3.txt', 'w') );

#subg_from_file = pickle.load(open('sub3.txt'));
#print "Form file\n", sub3;

genenameidx = 2; # default value = 0 (set it to 0 if user doesn't input any val)
foldchangeidx = 9; # default value = 1 # IF LOG FOLD CHANGE, ELIMINATE INF AND -INF
state1idx = 7; # NI
state2idx = 8; # ND # ALLOW USER TO INPUT DATA AS EXPRESSION VALUE OR RATIO

exp_obj = geneexp(expf, genenameidx, foldchangeidx, state1idx, state2idx); 
genelabels = exp_obj.assignNodeLabels(); # assign labels to gene based on gene expression
# print genelabels

pgenes = {}; # dict of dict: key: pathway name, value: dictionary where key is pathway gene id (entry id= node name) and value is gene label
allpathgenes = []; # list of all genes in all pathways: each entry is a list that has gene name and its aliases
#exit()

nodelabels = ['nn','n','p','pp'];
edgelabels = ['activate', 'inhibit'];
gdb = graphs();

def expLabelsPathgenes( pgenedict ):
    pgenelabels = {}; # this dict will be used to relabel pathway graph with gene-exp labels instead of entry-ids from kegg xml files
    print pgenedict

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
        print id, pgenenames, pglabel
        # if pathway genes have no label (i.e. no expression data) or multiple labels, remove those genes (nodes) from the pathway before generating subgraphs        
        if (len(pglabel) == 1): # each node/pgene should have exactly one label
            pgenelabels[id] = pglabel[0]; # assign gene exp label to entry id (node id)

    return pgenelabels;
	

# DO THIS ONLY THE FIRST TIME THE CODE IS RUN
# read all kgml files in KEGGpathways folder

num3graphs = 0; # number of trigraphs in all pathways
allkgmlfiles = [f for f in listdir("./KEGGpathways/")]

flag = 0
for f in allkgmlfiles:
    if ("sub" in f): # subgraphs have already been generated for the pathways
        flag = 1;
        print "file is ", f, " flag is ", flag;

if (flag == 0): # pathway files have not been read
    for f in allkgmlfiles:
        if (".xml" not in f):
            continue  # only read .xml files
        pobj = pathway("./KEGGpathways/"+f);
        pgenes[f] = expLabelsPathgenes( pobj.genedict );
        print "mapping: ", pgenes[f];
        print "Nodes before: ", pobj.pathwaygraph.nodes();
        pobj.removeNodes( pgenes[f] ); # remove all nodes that do not have a mapping to a gene exp label
        print "Nodes after: ", pobj.pathwaygraph.nodes();
        pobj.addExpLabelNode( pgenes[f] ); # add gene exp label as node attribute
        print NX.get_node_attributes(pobj.pathwaygraph, 'elabel')
        print pobj.pathwaygraph.edges(data=True);
        print pobj.pathwaygraph.nodes(data=True);
        #pathgenes =  pobj.genedict.values(); # all genes in this pathway
        sub3 = pobj.getSubgraphs( 3 ); # generate all connected subgraphs of size 3
        print sub3[0].edges(data=True);
        print sub3[0].nodes(data=True);
        pickle.dump(sub3, open("./KEGGpathways/"+f+'.sub3', 'w') );
        num3graphs += len(sub3);
        print f, len(sub3);
        gdb.populateGDB( sub3 );
        #gdb.getStats();
        # assign labels to pathway genes

print "Total number of trigraphs", num3graphs;

