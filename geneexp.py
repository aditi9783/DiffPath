''' Class to read gene expression data.'''

import numpy
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import pylab
import math
import collections

# change name to GeneExp 
class geneexp: 
    # rename to read_GeneExp
    def readGeneExp(self, fname, idx_gene, ifc, i2, i3):
        '''write a function to deal with replicate data'''
        #with open(fname, 'r') as fh:
	fh = open(fname, 'r')
        
	# remove first line from the file that has column headers: 
	# TEST IF THIS IS NEEDED (content[i2] should be float)
	fh.readline() 
        
	for line in fh:
           content = line.split('\t')
           after = float(content[i3])
           before = float(content[i2])
 	   # both genes don't show any expression: skip these because 
	   # a)not interesting, b) could be due to insufficient data
           if before == after and before == 0.0:
               continue;
           elif before == 0.0: # gene is turned on in the 'after' sample
               self.on_or_off_genes[ content[idx_gene] ] = "ON"
           elif after == 0.0: # gene is turned off in the 'after' sample
               self.on_or_off_genes[ content[idx_gene] ] = "OFF"
           else:
               self.geneexpdict[content[idx_gene]] = float(content[ifc]); # simple ratio of gene exp (after/before) is not normal, that's why use log fold change as a measure of expression change



    def statsAnalysis( self ):

    	def make_QQ( val_list, type_ ):
	    ge_arr = numpy.array( val_list )
            descstats = stats.describe( ge_arr ) # descriptive statistics for the log fold change values: size of array, (min,max), mean, var, skewness, kurtosis
            print descstats
            avg_logfc = numpy.mean( ge_arr );
            stdev_logfc = numpy.std( ge_arr );
            print type_, " mean and sd: ", avg_logfc, stdev_logfc;
            stats.probplot( ge_arr, plot=matplotlib.pyplot )
            matplotlib.pyplot.savefig('qqplot_'+type_+'.png')
            matplotlib.pyplot.close();
            return avg_logfc, stdev_logfc
 
	avg_logfc, stdev_logfc = make_QQ( self.geneexpdict.values(), "raw" )

        # if the distribution is not central, the n and nn labels could be assigned to genes with > 0 log(fc). To avoid this, convert gene exp values to z-scores and recalculate mean and sd
        for k, v in self.geneexpdict.items():
            zscore = (v - avg_logfc)/float(stdev_logfc);
            self.geneexpdict[k] = zscore;

        # recompute distribution parameters
	make_QQ( self.geneexpdict.values(), "centralized" )


    def assignNodeLabels( self ): # assign labels as: 'n' for negative (underexpressed), 'p' for positive (overexpressed)
    # -1 sd < log(fc) < mu => label 'n' for genes that are weakly underexpressed, mu and sd are from log(fc) distribution
    # log(fc) < -1 sd => label 'nn' for genes that are strongly underexpressed (below -1 SD)
    # mu < log(fc) < 1 sd => label 'p' for genes that are weakly overexpressed
    # log(fc) > 1 sd => label 'pp' for genes that are strongly overexpressed (above 1 SD)
    # these definitions can be further extended to 'nnn' for log(fc) < -2 sd, and 'ppp' for log(fc) > 2 sd              
    # Also leave genes whose expression is too close to the mean. Set lowerthres for that.

	def label_OnOffGenes():
            # genes that are turned off get 'nn' and those that are turned on get 'pp' label
            # not 'nnn' and 'ppp' because the on/off signal could just be due to less data
            for gene in self.on_or_off_genes:
                if (self.on_or_off_genes[gene] == "ON"):
                    genelabels[gene] = 'pp';
                    self.labeldist['pp'] += 1;
                elif (self.on_or_off_genes[gene] == "OFF"):
                    genelabels[gene] = 'nn';
                    self.labeldist['nn'] += 1;

	def assign_ExpLabel(log_fc, avg_change):
            if log_fc < avg_change:
	        label = 'n'
            elif log_fc > avg_change:
                label = 'p'
            else:
                return ""  

	    for i, thres in enumerate(thresholds):
                if (log_fc < avg_change-thres) or (log_fc > avg_change+thres):
                    return label*(i+1)
            
        genelabels = {};
        sigfac1 = 1 # how many SD away from mean should we assign labels pp and nn 
        sigfac2 = 2 # how many SD away from mean should we assign labels ppp and nnn 
        stdev_thres1 = sigfac1 * self.stdev_logfoldchange
        stdev_thres2 = sigfac2 * self.stdev_logfoldchange
        lowerthres = 0.01 * self.stdev_logfoldchange; # all genes with logfc lower than mean+lowerthres and higher than mean-lowerthres will be unlabeled and thus ignored (genes too close to average log fc)
        print "lowerthres: ", lowerthres, " stdevthres1: ", stdev_thres1, " stdevthres2: ", stdev_thres2;

	thresholds = [lowerthres, stdev_thres1, stdev_thres2]

        for gene, log_fc in self.geneexpdict.items():
            exp_label = assign_ExpLabel(log_fc, self.avg_logfoldchange)
            genelabels[gene] = exp_label
            self.labeldist[exp_label] += 1;
        print "Number of unlabeled genes: ", self.labeldist[""];

        label_OnOffGenes()
        return genelabels;

    def __init__ (self, expfile, genenameidx, foldchangeidx, beforeidx, afteridx):
        self.geneexpdict = {}; # dictionary with key = gene name, value = ratio of gene exp/ log fold change in gene exp. This does not include genes that are turned on or off in the 'after' sample. Those are stored separately.
        self.on_or_off_genes = {};
        self.avg_logfoldchange = 0.0; # initialize mean and sd for log fold change distribution
        self.stdev_logfoldchange = 0.0;
        self.labeldist = collections.Counter()

        self.readGeneExp( expfile, genenameidx, foldchangeidx, beforeidx, afteridx ) # index for the gene name and log2 fold change or ratio of expression 
        self.statsAnalysis(); # draws q-q plot for the distribution and computes, mean, sd for the gene-exp values
#        self.assignNodeLabels(); # label genes based on the gene expression
        print "Number of on/off genes: ", len(self.on_or_off_genes);
        print "Number of genes with expression value in both conditions: ", len(self.geneexpdict);
        #return self.genelabels; 
         
