PhyloBlend
==========

###**Description:**

Python pipeline to pull out sequence sets for the cluster-based mid-point rooted tree set released from the Phylota database. It will also allow expansion of the sequence set with BLAST searches, followed up with a new sets of alignments and tree construction. Other features are in there, but not described here yet.

<i>**Note:** Instructions below are for a [Debian](www.debian.org) based OS like [Ubuntu](www.ubuntu.com)</i>

#####**Requirements (Quick List):**


1. [Python 2.7](www.python.org)
2. [Numpy](www.numpy.org)
3. [Sci-py](www.scipy.org)
4. [Biopython](www.biopython.org/wiki/Main_Page)
5. [IVY](www.reelab.net/home/software/ivy/)
6. [Muscle](www.drive5.com/muscle/) for Alignments
7. [RAxML](sco.h-its.org/exelixis/web/software/raxml/index.html) for Tree Building
8. Local copy of [Phylota](www.phylota.net/pb/Download/) mySQL database
9. Local copy of [NCBI BLAST](www.blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download) 
10. Local copy of the NCBI [*nt*](ftp.ncbi.nlm.nih.gov/blast/db) database

***

####**Detailed Setup:**

####Install the necessary requirements:

1. Install [Numpy](www.numpy.org) <pre>sudo apt-get install python-numpy</pre>
2. Install [Sci-py](www.scipy.org) <pre>sudo apt-get install python-scipy</pre>
3. Install [Biopython](www.biopython.org/wiki/Main_Page) <pre>sudo apt-get install python-biopython</pre>
4. Install [IVY](www.reelab.net/home/software/ivy/)
````
1. git clone https://github.com/rhr/ivy
2. sudo python setup.py install
````
5. Install [Muscle](www.drive5.com/muscle/) for Alignments <pre>sudo apt-get install muscle</pre>
6. Install [RAxML](sco.h-its.org/exelixis/web/software/raxml/index.html) for Tree Building <pre>sudo apt-get install raxml</pre>
7. Download local copy of the [Phylota](www.phylota.net/pb/Download/) mySQL database dump <pre> code here</pre>
8. Prime MySQL with the Phylota data <pre>commands here</pre>
9. Download a local copy of [NCBI BLAST](www.blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)<pre>sudo apt-get install ncbi-blast+</pre>  
10. Download a local copy of the NCBI [*nt*](ftp.ncbi.nlm.nih.gov/blast/db) database <pre> code here </pre>

***

####**Example Use:**

Work in Progress


