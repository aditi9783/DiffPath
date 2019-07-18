# DiffPath

DiffPath is a program to detect differentially expressed sub-pathways from RNAseq and pathway (in xml format, as in KEGG database) data.

Benefits of using DiffPath:
- Takes pathway structure into account (unlike GSEA)
- Can detect a broader set of differentially expressed sub-pathways (unlike SPIA that is limited to tree-like structure sub-pathways)
- Is not biased towards most significantly differentially expressed genes in a pathway (unlike DEAP)

DiffPath works in three steps (see visual in diffpath_alg.pptx):
- Binning raw gene-expression values, thus discretizing the continuous data into a finite set of labels. This reduces bias towards genes with extreme gene-expression values.
- Generating all sub-pathways of desired size (size=number of genes, or nodes, in sub-pathway) and applying gene-expression labels to the genes in the sub-pathways. The pathway structure (different edge-types for activating and inhibitory gene interactions) is left intact. 
- Comparing each labeled sub-pathway to a thousand randomizations of node-labels in the entire set of sub-pathways and reporting the ones that are significantly differentially expressed at the desired FDR level.
