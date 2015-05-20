# functions for graph enumeration and graph matching

import networkx as NX
import itertools

class graphs:
    """
    This class enumerates and indexes graphs, does graph matching, and generate statistics 
    """
    def getNodeEdgeKeys( self, ndict, edict ): # create keys like keys in connectivitydict, for all edges in subgraph
        keys = []; # list of all node-edge keys in the graph
        for nodes,edge in edict.iteritems():
            keys.append( "-".join( [ndict[nodes[0]], ndict[nodes[1]], edge] ) );
        return keys;   

    def getGraphID( self, keys, connectivitydict ): # generate an index of length 30, mapping frequncies of node-edge keys in the subgraph
        idx = [0] * 63; # initialize graph id
        for k in keys:
            kid = connectivitydict[k];
            idx[ kid ] += 1 # increase count of node-edge key in graph id
        return "".join(str(x) for x in idx); # return string rep of graph

    def populateGDB( self, glist, nlabels, connectivitydict ): # take in a list of subgraphs and populate the gdb with it
        gidmap = {}; # key is gid (graph id from nodes and edges), and value is list of subgraphs with that gid
#        print "number of graphs: ", len(glist);
        for g in glist:
            # get node labels as dict that has key: node id, value: expression label (p,n,pp,nn)
            # nlabels = NX.get_node_attributes(g, 'elabel');
            # get edge labels as dict that has key: ids of nodes that have the edge, value: edge type
            edgelabels = NX.get_edge_attributes(g, 'type');
            #print nlabels;
            #print edgelabels;
            gkeys = self.getNodeEdgeKeys(nlabels, edgelabels);
            gid = self.getGraphID( gkeys, connectivitydict );
#            print gid
            if gid in self.gdb: # if gid for this subgraph already exists in gdb, increment its count, else create new gid
                self.gdb[gid] += 1;
                gidmap[gid].append( g );
            else:
                self.gdb[gid] = 1;
                gidmap[gid] = [g]; 
        return gidmap;

    def getStats( self ):
        num_gids = len( self.gdb.values() );
        print "Number of distinct subgraphs: ", num_gids;
        # print self.gdb.values();
        print "Total number of subgraphs: ", sum(self.gdb.values());

    def __init__( self ): # read kegg xml file and parse it
        """
        Parameters: nothing

        Returns: nothing

        Description: Initializes graph data base (gdb): indexed graphs for graph matching and computing statistics 
        """
        self.gdb = {};
        #self.initializeGDB(); # keys: binary rep of node lables in subgraphs
        # 'p','n','pp','nn'. Thus subgraphs with three 'p' node lables will have key 1000
