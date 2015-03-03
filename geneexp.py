# Class to read gene expression data

class geneexp: 

    def readGeneExp( self, fname, i1, i2 ): # write a function to deal with replicate data
        fh = open(fname, 'r')
        fh.readline()  # remove first line from the file that has column headers
        for line in fh:
           content = line.split('\s+')
           if (float(content[i1]) > 0.0 or float(content[i2]) > 0.0):
               yield content # return line only if at least one of the genes is showing some expression


    def __init__ (self, expfile, idx1, idx2):
        self.validexp = self.readGeneExp( expfile, idx1, idx2 ) # indices for the two conditions: first one is default and second is perturbed
        print self.validexp
