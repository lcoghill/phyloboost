Phyloboost
==========

###**Description:**

Python pipeline to pull out sequence sets for the cluster-based mid-point rooted tree set released from the Phylota database. It can handle the entire tree set, or any subset of trees. It will also allow expansion of the sequence set with BLAST searches, followed up with a new sets of alignments and tree construction. Other features are coming soon.

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
11. [Git]()

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
9. Install and Prime MySQL with the Phylota data<pre>
$ sudo apt-get install mysql-server
$ mysql -uUSERHERE -p  --max_allowed_packet=300M --connect_timeout=6000
mysql> CREATE DATABASE phylota;
mysql> source pb.bu.rel184.4.10.2012.gz
mysql> quit
</pre>
*Note:* This will take quite awhile.
10. Download the Genbank nt nucleotide set<pre>
$ wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nt.gz</pre>
11. Convert the nt fasta file to csv<pre>
$ python fasta-csv.py</pre>
12. Import the nt csv file into mysql<pre>
$ mysql -uUSERHERE -p  --max_allowed_packet=500M --connect_timeout=6000 --local-infile=1
mysql> SELECT DBNAMEHERE;
mysql> DROP TABLE IF EXISTS genbank_nucleotide;
mysql> CREATE TABLE genbank_nucleotide(id INT PRIMARY KEY AUTO_INCREMENT, gi INT(20), accession VARCHAR(25), length int(10), title VARCHAR(25), sequence LONGTEXT) ENGINE=MyISAM DEFAULT CHARSET=latin1 AUTO_INCREMENT=1;
mysql> LOAD DATA LOCAL INFILE '/full/path/here/nt.fasta.csv' INTO TABLE genbank_nucleotide FIELDS TERMINATED BY ',' lines terminated by '\n' (gi, accession, length, title, sequence);
mysql> quit
</pre>
*Note:* This will take quite awhile.
13. Download a local copy of [NCBI BLAST](http://www.blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)<pre>sudo apt-get install ncbi-blast+</pre>  
14. Download a local copy of the NCBI [*nt*](ftp://ftp.ncbi.nlm.nih.gov/blast/db) database <pre>
$ wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.2.29+-x64-linux.tar.gz
$ tar -zxvf ncbi-blast-2.2.29+-x64-linux.tar.gz
$ mkdir ~/blastdb
$ cp ncbi-blast-2.2.29+-x64-linux.tar.gz/bin/update_blastdb.pl ~/blastdb
$ perl update_blastdb.pl --passive --decompress nt
</pre>
*Note:* This will take awhile due to the size and number of files.
15. Install Git <pre>sudo apt-get install git</pre>
16. Clone this repository into desired location <pre>git clone https://github.com/lcoghill/phyloboost</pre>

***

####**Example Use:**

Work in Progress


