# functions for graph enumeration and graph matching

import networkx as NX
import itertools

nodelabeldict = {'p' : 0, 'n' : 1, 'pp' : 2, 'nn' : 3}; # index of node label in gdb key
edgelabeldict = {'activate' : 0, 'inhibit' : 1, 'interact' : 2, 'none' : 3}; # index of edge label in gdb
connectivitydict = { 'p-p-activate' : 0,
                     'p-n-activate' : 1,
                     'n-p-activate' : 1,
                     'p-pp-activate' : 2,
                     'pp-p-activate' : 2,
                     'p-nn-activate' : 3,
                     'nn-p-activate' : 3,
                     'n-n-activate' : 4,
                     'n-pp-activate' : 5,
                     'pp-n-activate' : 5,
                     'n-nn-activate' : 6,
                     'nn-n-activate' : 6,
                     'pp-pp-activate' : 7,
		     'pp-nn-activate' : 8,
		     'nn-pp-activate' : 8,
                     'nn-nn-activate' : 9,
                     'p-p-inhibit' : 10,
                     'p-n-inhibit' : 11,
                     'n-p-inhibit' : 11,
                     'p-pp-inhibit' : 12,
                     'pp-p-inhibit' : 12,
                     'p-nn-inhibit' : 13,
                     'nn-p-inhibit' : 13,
                     'n-n-inhibit' : 14,
                     'n-pp-inhibit' : 15,
                     'pp-n-inhibit' : 15,
                     'n-nn-inhibit' : 16,
                     'nn-n-inhibit' : 16,
                     'pp-pp-inhibit' : 17,
		     'pp-nn-inhibit' : 18,
		     'nn-pp-inhibit' : 18,
                     'nn-nn-inhibit' : 19,
                     'p-p-interact' : 20,
                     'p-n-interact' : 21,
                     'n-p-interact' : 21,
                     'p-pp-interact' : 22,
                     'pp-p-interact' : 22,
                     'p-nn-interact' : 23,
                     'nn-p-interact' : 23,
                     'n-n-interact' : 24,
                     'n-pp-interact' : 25,
                     'pp-n-interact' : 25,
                     'n-nn-interact' : 26,
                     'nn-n-interact' : 26,
                     'pp-pp-interact' : 27,
		     'pp-nn-interact' : 28,
		     'nn-pp-interact' : 28,
                     'nn-nn-interact' : 29 };

class graphs:
    """
    This class enumerates and indexes graphs, does graph matching, and generate statistics 
    """
    def getLabelKey( self, list, dict ):
        lkey = ["0","0","0","0"];
        for l in list:
            lidx = dict[l];
            lkey[lidx] = "1";
        return "".join(lkey);

    def getNodeEdgeKeys( self, ndict, edict ): # create keys like keys in connectivitydict, for all edges in trigraph
        keys = []; # list of all node-edge keys in the graph
        for nodes,edge in edict.iteritems():
            keys.append( "-".join( [ndict[nodes[0]], ndict[nodes[1]], edge] ) );
        return keys;   

    def getGraphID( self, keys ): # generate an index of length 30, mapping frequncies of node-edge keys in the trigraph
        idx = [0] * 30; # initialize graph id
        for k in keys:
            kid = connectivitydict[k];
            idx[ kid ] += 1 # increase count of node-edge key in graph id
        return "".join(str(x) for x in idx); # return string rep of graph

    def populateGDB( self, glist, nlabels ): # take in a list of trigraphs and populate the gdb with it
        gidmap = {}; # key is gid (graph id from nodes and edges), and value is list of trigraphs with that gid
#        print "number of graphs: ", len(glist);
        for g in glist:
            # get node labels as dict that has key: node id, value: expression label (p,n,pp,nn)
            # nlabels = NX.get_node_attributes(g, 'elabel');
            # get edge labels as dict that has key: ids of nodes that have the edge, value: edge type
            edgelabels = NX.get_edge_attributes(g, 'type');
            #print nlabels;
            #print edgelabels;
            gkeys = self.getNodeEdgeKeys(nlabels, edgelabels);
            gid = self.getGraphID( gkeys );
#            print gid
            if gid in self.gdb: # if gid for this trigraph already exists in gdb, increment its count, else create new gid
                self.gdb[gid] += 1;
                gidmap[gid].append( g );
            else:
                self.gdb[gid] = 1;
                gidmap[gid] = [g]; 
        return gidmap;

    def getStats( self ):
        num_gids = len( self.gdb.values() );
        print "Number of distinct trigraphs: ", num_gids;
        print self.gdb.values();
        print "Total number of subgraphs: ", sum(self.gdb.values());

    def __init__( self ): # read kegg xml file and parse it
        """
        Parameters: nothing

        Returns: nothing

        Description: Initializes graph data base (gdb): indexed graphs for graph matching and computing statistics 
        """
        self.gdb = {};
        #self.initializeGDB(); # keys: binary rep of node lables in trigraphs
        # 'p','n','pp','nn'. Thus trigraphs with three 'p' node lables will have key 1000
