#!/bin/bash
## Downloads the genbank nucleotide dataset in fasta format
## and builds a local BLAST Database of eukaryotes for Phyloboost

$title = "Phyloboost version 1.5"
$version = "phyloboost_1.5"



wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nt.gz
gunzip nt.gz
mv nt nt.fas
makeblastdb -in nt.fas -out $version -dbtype 'nucl' -title $title -parse_seqids
mkdir $version
mv *.nsq $version
mv *.nhr $version
mv *.nin $version
mv *.nal $version