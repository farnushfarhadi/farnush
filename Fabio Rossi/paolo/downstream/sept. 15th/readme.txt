
Receptor 1 -> gene1 is the only think I reported in my database.

Concerning your request to have chains of genes, I am developing a way to implement it.
All the information needed are stored in graphite R package.

The thing I miss is the end point, let's say gene5. Unless we define the end point it is pretty impossible (for me at least :) ) to get all the shortest paths between receptors and all genes.

My question is: do I have to set as end point the transcription factors? what else?

I made some effort (please do not consider them complete) in this directions:

a) I downloaded a list of TFs for Mus musculus.

b) I checked those expressed exclusively in one single cell (just like the analysis for ligands and receptor)

c) I build a meta-netowork with all the interaction (edges in kegg) that refer to "Process(activation)" or "Process(expression)"; I made the network like this to keep the signal as straight as possible. In other word any arrow in the chain gene1->gene2->gene3... so on represents either a "Process(activation)" or "Process(expression)".

d) I took the receptor from the previous analysis (the one I sent one week ago more or less)

e) I computed all the shortest paths between all receptor to TF pairs.

The results is very big also considering that the TF number was dramatically reduced! This results make me consider the possibility to do this think only on few pathways (even if it is pretty long) and on single pathway.

In attachment I send you 3 file: one for every cell-type.

The 3 files have exactly the same format:
column 1 -> receptor symbol
column 2 -> tf symbol
column 3 -> chain of signal from the 1st interactor of the receptor to the tf itself