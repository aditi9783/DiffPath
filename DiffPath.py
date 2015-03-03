
import sys
from pathway import *
from geneexp import *
import xml.etree.ElementTree as ET

pathwayf = sys.argv[1];
expf = sys.argv[2];

#path_obj = pathway(pathwayf); # read kgml pathway file
#obj.getSubPathways(); # get connected components in the pathway file

geneexp(expf, 1, 2);
