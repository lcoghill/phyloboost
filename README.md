PhyloBlend
==========

###**Description:**

Python pipeline to pull out sequence sets for the cluster-based mid-point rooted tree set released from the Phylota database. It will also allow expansion of the sequence set with BLAST searches, followed up with a new sets of alignments and tree construction. Other features are in there, but not described here yet.

<i>**Note:** Instructions below are for a [Debian](http://www.debian.org) based OS like [Ubuntu](http://www.ubuntu.com)</i>

#####**Requirements (Quick List):**


1. [Python 2.7](http://www.python.org)
2. [Numpy](http://www.numpy.org)
3. [Sci-py](http://www.scipy.org)
4. [Biopython](http://www.biopython.org/wiki/Main_Page)
5. [IVY](http://www.reelab.net/home/software/ivy/)
6. [Muscle](http://www.drive5.com/muscle/) for Alignments
7. [RAxML](http://sco.h-its.org/exelixis/web/software/raxml/index.html) for Tree Building
8. Local copy of [Phylota](http://www.phylota.net/pb/Download/) mySQL database
9. Local copy of [NCBI BLAST](http://www.blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download) 
10. Local copy of the NCBI [*nt*](ftp.ncbi.nlm.nih.gov/blast/db) database

***

####**Detailed Setup:**

####Install the necessary requirements:

1. Install [Numpy](http://www.numpy.org) <pre>$ sudo apt-get install python-numpy</pre>
2. Install [Sci-py](http://www.scipy.org) <pre>$ sudo apt-get install python-scipy</pre>
3. Install [Biopython](http://www.biopython.org/wiki/Main_Page) <pre>$ sudo apt-get install python-biopython</pre>
4. Install [IVY](http://www.reelab.net/home/software/ivy/)<pre>
$ git clone https://github.com/rhr/ivy
$ sudo python setup.py install
</pre>
5. Install [Muscle](http://www.drive5.com/muscle/) for Alignments <pre>$ sudo apt-get install muscle</pre>
6. Install [RAxML](http://sco.h-its.org/exelixis/web/software/raxml/index.html) for Tree Building <pre>$ sudo apt-get install raxml</pre>
7. Download local copy of the [Phylota](http://www.phylota.net/pb/Download/) mySQL database dump <pre>$ wget http://www.phylota.net/pb/Download/184/pb.bu.rel184.4.10.2012.partaa
$ wget http://www.phylota.net/pb/Download/184/pb.bu.rel184.4.10.2012.partab
$ wget ...partae
</pre>
8. Combine the pieces into a single .sql file<pre>$ cat pb.bu.rel184.4.10.2012.part* > pb.bu.rel184.4.10.2012.gz</pre>
9. Prime MySQL with the Phylota data <pre>commands here</pre>
10. Download a local copy of [NCBI BLAST](http://www.blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)<pre>sudo apt-get install ncbi-blast+</pre>  
11. Download a local copy of the NCBI [*nt*](ftp://ftp.ncbi.nlm.nih.gov/blast/db) database <pre>
$ wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.2.29+-x64-linux.tar.gz
$ tar -zxvf ncbi-blast-2.2.29+-x64-linux.tar.gz
$ mkdir ~/blastdb
$ cp ncbi-blast-2.2.29+-x64-linux.tar.gz/bin/update_blastdb.pl ~/blastdb
$ perl blast_updatedb.pl --passive --decompress nt
</pre>
*Note:* This will take awhile due to the size and number of files.
12. 
***

####**Example Use:**

Work in Progress


