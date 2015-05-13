# Class to read gene expression data

import numpy
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import pylab
import math

class geneexp: 

    def readGeneExp( self, fname, i1, ifc, i2, i3): # write a function to deal with replicate data
        fh = open(fname, 'r')
        fh.readline()  # remove first line from the file that has column headers: TEST IF THIS IS NEEDED (content[i2] should be float)
        for line in fh:
           content = line.split('\t')
           #print content
           after = float( content[i3] )
           before = float( content[i2] )
           if ( before == after and before == 0.0 ): # both genes don't show any expression: skip these because a)not interesting, b) could be due to insufficient data
               continue;
           elif ( before == 0.0 ): # gene is turned on in the 'after' sample
               self.on_or_off_genes[ content[i1] ] = "ON"
           elif ( after == 0.0 ): # gene is turned off in the 'after' sample
               self.on_or_off_genes[ content[i1] ] = "OFF"
           else:
               self.geneexpdict[content[i1]] = float(content[ifc]); # simple ratio of gene exp (after/before) is not normal, that's why use log fold change as a measure of expression change
               #self.geneexpdict[content[i1]] = after/before; # ratio of gene exp: after/before

    def statsAnalysis( self ):
        ge_arr = numpy.array( self.geneexpdict.values() )
        descstats = stats.describe( ge_arr ) # descriptive statistics for the log fold change values: size of array, (min,max), mean, var, skewness, kurtosis
        print descstats
        self.avg_logfoldchange = numpy.mean( ge_arr );
        self.stdev_logfoldchange = numpy.std( ge_arr );
        print "mean and sd: ", self.avg_logfoldchange, self.stdev_logfoldchange;
        stats.probplot( ge_arr, plot=matplotlib.pyplot )
        matplotlib.pyplot.savefig('qqplot.png')

    def assignNodeLabels( self ): # assign labels as: 'n' for negative (underexpressed), 'p' for positive (overexpressed)
    # -1 sd < log(fc) < mu => label 'n' for genes that are weakly underexpressed, mu and sd are from log(fc) distribution
    # log(fc) < -1 sd => label 'nn' for genes that are strongly underexpressed (below -1 SD)
    # mu < log(fc) < 1 sd => label 'p' for genes that are weakly overexpressed
    # log(fc) > 1 sd => label 'pp' for genes that are strongly overexpressed (above 1 SD)
    # these definitions can be further extended to 'nnn' for log(fc) < -2 sd, and 'ppp' for log(fc) > 2 sd              
        genelabels = {};
        sigfac = 1 # how many SD away from mean should we assign labels for strong change in exp
        stdev_thres = sigfac * self.stdev_logfoldchange
        for gene in self.geneexpdict:
            log_fc = self.geneexpdict[gene];
            if (log_fc < self.avg_logfoldchange):
                if (log_fc < self.avg_logfoldchange-stdev_thres):
                    genelabels[gene] = 'nn'
                else:
                    genelabels[gene] = 'n'     
            elif (log_fc > self.avg_logfoldchange):
                if (log_fc > self.avg_logfoldchange+stdev_thres):
                    genelabels[gene] = 'pp'
                else:
                    genelabels[gene] = 'p'
            # if gene expression log(fc) is exactly at mean, ignore that gene
            # print gene, log_fc, genelabels[gene];
        # genes that are turned off get 'nn' and those that are turned on get 'pp' label
        for gene in self.on_or_off_genes:
            if (self.on_or_off_genes[gene] == "ON"):
                 genelabels[gene] = 'pp'
            elif (self.on_or_off_genes[gene] == "OFF"):
                 genelabels[gene] = 'nn'
        return genelabels

    def __init__ (self, expfile, genenameidx, foldchangeidx, beforeidx, afteridx):
        self.geneexpdict = {}; # dictionary with key = gene name, value = ratio of gene exp/ log fold change in gene exp. This does not include genes that are turned on or off in the 'after' sample. Those are stored separately.
        self.on_or_off_genes = {};
        self.avg_logfoldchange = 0.0; # initialize mean and sd for log fold change distribution
        self.stdev_logfoldchange = 0.0;

        self.readGeneExp( expfile, genenameidx, foldchangeidx, beforeidx, afteridx ) # index for the gene name and log2 fold change or ratio of expression 
        self.statsAnalysis(); # draws q-q plot for the distribution and computes, mean, sd for the gene-exp values
#        self.assignNodeLabels(); # label genes based on the gene expression
        print len(self.on_or_off_genes);
        #return self.genelabels; 
         
