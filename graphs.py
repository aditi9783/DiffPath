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

    def getGraphID( self, ndict, edict ):
        gid = []; # list of all node-edge keys in the graph
        for nodes,edge in edict.iteritems():
            gid.append( "-".join( [ndict[nodes[0]], ndict[nodes[1]], edge] ) );
        return gid;   

    def populateGDB( self, glist ): # take in a list of trigraphs and populate the gdb with it
        print "number of graphs: ", len(glist);
        for g in glist:
            # get node labels as dict that has key: node id, value: expression label (p,n,pp,nn)
            nlabels = NX.get_node_attributes(g, 'elabel');
            # get edge labels as dict that has key: ids of nodes that have the edge, value: edge type
            edgelabels = NX.get_edge_attributes(g, 'type');
            print nlabels;
            print edgelabels;
            gid = self.getGraphID(nlabels, edgelabels);
            print gid;

    def populateGDB_0( self, glist ): # take in a list of trigraphs and populate the gdb with it
        print "number of graphs: ", len(glist);
        for g in glist:
            # get node labels as values of the dict that has key: node id, value: expression label (p,n,pp,nn)
            nlabels = NX.get_node_attributes(g, 'elabel').values();
            # get edge labels as values of the dict that has key: ids of nodes that have the edge, value: edge type
            edgelabels = NX.get_edge_attributes(g, 'type').values();
            nodekey = self.getLabelKey( nlabels, nodelabeldict );
            edgekey = self.getLabelKey( edgelabels, edgelabeldict );
            print nlabels, nodekey;
            print edgelabels, edgekey;
            # add graph to the graph database
            self.gdb[nodekey][edgekey].append( g );
            print "added to ", nodekey+edgekey;
            print "\n", len(self.gdb["1100"]["1010"]), "\n";
        # print self.gdb;

    def getStats( self ):
        graphdist = {};
        for nkey in self.gdb:
            for ekey in self.gdb[nkey]:
                graphdist[nkey+ekey] = len(self.gdb[nkey][ekey]); # number of graphs with that combination of node and edge keys
        for k,v in graphdist.iteritems():
            print k,v;

    def initializeGDB( self ):
    # Two layered hash: first key represents node labels in each trigraph, second key represents edgelabels, the value of second key holds list of trigraphs with the edge and node lables given by the first two keys
        gdbkeys = ["".join(seq) for seq in itertools.product("01", repeat=4)]; # generates all binary seq of length 4 (all possible instances of 4 node labels)
        edgekeys = gdbkeys;
        edgedict = dict( zip(edgekeys, [ [] for i in range(0,len(edgekeys)) ] ) ); # keys: binary seq for edges in trigraph, val: list of all trigraphs that match that edge set
        # remove edgekey "0001" (there won't be a graph with all 'none' edges); "1111" (no trigraphs with all three edge types and also 'none' edgetype); 0000 (no trigraphs with no edges)
        edgedict.pop("0001",None);
        edgedict.pop("1111",None);
        edgedict.pop("0000",None);

        gdbvals = [ edgedict for i in range(0,len(gdbkeys)) ];
        self.gdb = dict( zip(gdbkeys, gdbvals) );
        # remove the '0000' key as at least one node label is expected to be there (even if all three nodes have same label)
        # remove the '1111' key as trigraphs can have max 3 nodel labels 
        self.gdb.pop("0000",None);
        self.gdb.pop("1111",None);

    def initializeGDB_0( self ):
    # Two layered hash: first key represents node labels in each trigraph, second key represents edgelabels, the value of second key holds list of trigraphs with the edge and node lables given by the first two keys
        gdbkeys = ["".join(seq) for seq in itertools.product("01", repeat=4)]; # generates all binary seq of length 4 (all possible instances of 4 node labels)
        edgekeys = gdbkeys;
        edgedict = dict( zip(edgekeys, [ [] for i in range(0,len(edgekeys)) ] ) ); # keys: binary seq for edges in trigraph, val: list of all trigraphs that match that edge set
        # remove edgekey "0001" (there won't be a graph with all 'none' edges); "1111" (no trigraphs with all three edge types and also 'none' edgetype); 0000 (no trigraphs with no edges)
        edgedict.pop("0001",None);
        edgedict.pop("1111",None);
        edgedict.pop("0000",None);

        gdbvals = [ edgedict for i in range(0,len(gdbkeys)) ];
        self.gdb = dict( zip(gdbkeys, gdbvals) );
        # remove the '0000' key as at least one node label is expected to be there (even if all three nodes have same label)
        # remove the '1111' key as trigraphs can have max 3 nodel labels 
        self.gdb.pop("0000",None);
        self.gdb.pop("1111",None);

        # print self.gdb;

    def __init__( self ): # read kegg xml file and parse it
        """
        Parameters: nothing

        Returns: nothing

        Description: Initializes graph data base (gdb): indexed graphs for graph matching and computing statistics 
        """
        self.gdb = {};
        #self.initializeGDB(); # keys: binary rep of node lables in trigraphs
        # 'p','n','pp','nn'. Thus trigraphs with three 'p' node lables will have key 1000
