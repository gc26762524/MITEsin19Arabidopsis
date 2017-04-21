# MITEs in 19 Arabidopsis Accessions #
This repository contains the MySQL database used in "Genome-Wide Comparative Analysis of Miniature Inverted Repeat Transposable Elements in 19 *Arabidopsis thaliana* Ecotype Accessions." Both a zip file containing a dump of the database and the scripts to populate the database from scratch are available in this repository. This tutorial is for running on linux (Ubuntu).

## Installing MySQL ##

To install MySQL locally run the following:

```
sudo apt-get install mysql-server
```

In order to run the scripts to generate the database, in the source code and this tutorial the username and root have been set to an example of 'root' and 'BioInfoLab' respectively. Using a real username and password will require updating the source code. Specifically in mitedb.py look for the credential string:

```
mysql+pymysql://root:BioInfoLab@localhost
```

Additionally, if using a MySQL server running on a different machine (instead of locally), change `@localhost` in the above credential string to point to that server instead. However, this is not necessary when restoring the database from the provided zip file.

## Restoring the Database from the Zip File ##

First unzip the file:

```
unzip MITEdb.zip
```

Restore the database into mysql:

```
mysql -u root -p
mysql> create database MITEdb;
mysql> use MITEdb;
mysql> source MITEdb.sql;
```

## Generating the Database from Source ##

In this tutorial, version numbers of dependencies resolved via `conda` and `pip` commands will be specified to match the environment at the time the code was written and used.

### Python Environment Configuration ###

To get started quickly download and install the latestest copy [miniconda for Python 3.X](http://conda.pydata.org/miniconda.html) if miniconda or anaconda are not already on the system. Navigate to the `source` folder create a python environment:

```
conda create -n mitedb biopython=1.65 nose=1.3.7 pip=7.1.2 sqlalchemy=1.0.8 python=3.4
```

Active the environment:

```
source activate mitedb
```

Run these commands to finish configuring the mitedb Python environment:
   
```
pip install -Iv bcbio-gff==0.6.2
pip install -Iv PyMySQL==0.6.7
```

Skip this next step if there are no errors. If you get a setuptools error for any of the above try running the following command and then try the pip commands again:

```
conda install -f setuptools=18.3.2
```

### Populating the database using the Python Scripts ###

This project had a few unit tests, verify that they pass:

```
nosetests test.py
```

The following will run the main python script to create the database:

```
python mitedb.py
```

The ouput should show something like this:

```
Creating mysql db (overriding if exists)... Finished
Reading MITE/Gene correlations from bur.CalculationOfMITEsAndGenes... Finished
Reading MITE/Gene correlations from col.CalculationOfMITEsAndGenes... Finished
Reading MITE/Gene correlations from can.CalculationOfMITEsAndGenes... Finished
Reading MITE/Gene correlations from kn.CalculationOfMITEsAndGenes... Finished
Reading MITE/Gene correlations from edi.CalculationOfMITEsAndGenes... Finished
Reading MITE/Gene correlations from wu.CalculationOfMITEsAndGenes... Finished
Reading MITE/Gene correlations from hi.CalculationOfMITEsAndGenes... Finished
Reading MITE/Gene correlations from no.CalculationOfMITEsAndGenes... Finished
Reading MITE/Gene correlations from ws.CalculationOfMITEsAndGenes... Finished
Reading MITE/Gene correlations from rsch.CalculationOfMITEsAndGenes... Finished
Reading MITE/Gene correlations from zu.CalculationOfMITEsAndGenes... Finished
Reading MITE/Gene correlations from ct.CalculationOfMITEsAndGenes... Finished
Reading MITE/Gene correlations from oy.CalculationOfMITEsAndGenes... Finished
Reading MITE/Gene correlations from po.CalculationOfMITEsAndGenes... Finished
Reading MITE/Gene correlations from sf.CalculationOfMITEsAndGenes... Finished
Reading MITE/Gene correlations from wil.CalculationOfMITEsAndGenes... Finished
Reading MITE/Gene correlations from ler.CalculationOfMITEsAndGenes... Finished
Reading MITE/Gene correlations from tsu.CalculationOfMITEsAndGenes... Finished
Reading MITE/Gene correlations from mt.CalculationOfMITEsAndGenes... Finished
Importing Genes from denovo_annotation.rsch.gff3... Finished
Importing Genes from denovo_annotation.kn.gff3... Finished
Importing Genes from denovo_annotation.po.gff3... Finished
Importing Genes from denovo_annotation.zu.gff3... Finished
Importing Genes from denovo_annotation.tsu.gff3... Finished
Importing Genes from denovo_annotation.hi.gff3... Finished
Importing Genes from denovo_annotation.mt.gff3... Finished
Importing Genes from denovo_annotation.ct.gff3... Finished
Importing Genes from denovo_annotation.can.gff3... Finished
Importing Genes from denovo_annotation.sf.gff3... Finished
Importing Genes from denovo_annotation.wu.gff3... Finished
Importing Genes from denovo_annotation.col.gff3... Finished
Importing Genes from denovo_annotation.edi.gff3... Finished
Importing Genes from denovo_annotation.ler.gff3... Finished
Importing Genes from denovo_annotation.no.gff3... Finished
Importing Genes from denovo_annotation.oy.gff3... Finished
Importing Genes from denovo_annotation.wil.gff3... Finished
Importing Genes from denovo_annotation.ws.gff3... Finished
Importing Genes from denovo_annotation.bur.gff3... Finished
Importing MITEs from bur.fas.out.gff... Finished
Importing MITEs from col.fas.out.gff... Finished
Importing MITEs from can.fas.out.gff... Finished
Importing MITEs from kn.fas.out.gff... Finished
Importing MITEs from edi.fas.out.gff... Finished
Importing MITEs from wu.fas.out.gff... Finished
Importing MITEs from hi.fas.out.gff... Finished
Importing MITEs from no.fas.out.gff... Finished
Importing MITEs from ws.fas.out.gff... Finished
Importing MITEs from rsch.fas.out.gff... Finished
Importing MITEs from zu.fas.out.gff... Finished
Importing MITEs from ct.fas.out.gff... Finished
Importing MITEs from oy.fas.out.gff... Finished
Importing MITEs from po.fas.out.gff... Finished
Importing MITEs from sf.fas.out.gff... Finished
Importing MITEs from wil.fas.out.gff... Finished
Importing MITEs from ler.fas.out.gff... Finished
Importing MITEs from tsu.fas.out.gff... Finished
Importing MITEs from mt.fas.out.gff... Finished
```

When finished working in the environment run:

```
source deactivate mitedb
```

## Viewing the Database ##

Log in and select the MITEdb database:

```
mysql -u root -p
use MITEdb;
```

List the four tables:

```
show tables;
+------------------+
| Tables_in_MITEdb |
+------------------+
| gene             |
| gene_sub_feature |
| genome           |
| mite             |
+------------------+
```

The tables store the following information:

genome:
 - id (autogenerated)
 - accession (NCBI accession)
 - taxonomy (e.g. 3702)
 - code (e.g. 'bur', 'col', from file name of each of the 19 genome files)
 - comments

mite:
 - name (Assigned using the AT<chromosome>M<{10, 20, 30, ... n}> pattern. The same name in one accession does not necessarily represent the same MITE found in another accession even if they have the same name. This field along with genome_id are the primary key.)
 - chromosome
 - strand
 - start
 - end
 - correlated_gene (NULL if it is not known to be within 40 bps of an annotated or novel gene, essentially this affords us a mapping between mites and genes)
 - tir (the length of the TIR from the clustered reference MITE it belongs to)
 - tsd (the length of the TSD from the clustered reference MITE it belongs to)
 - family (from the clustered reference MITE it belongs to)
 - superfamily (from the clustered reference MITE it belongs to)
 - sequence (currently empty)
 - length
 - genome_id

gene:
 - name (Annotated name or "novel_gene" number from input files. This field along with genome_id are the primary keys)
 - chromosome
 - strand
 - start
 - end
 - exon_count (number of exons in the gene)
 - genome_id

gene_sub_features:
 - name (autogenerated, combined with genome_id serves as primary key)
 - gene_name (of the gene the sub feature belongs to)
 - type (one of: EXON, 5' UTR, 3' UTR)
 - start
 - end
 - genome_id (of the genome/accession the sub feature belongs to)

Example Query:

```
select mite.name, genome.code, mite.start,
       mite.end, sub.type, sub.start as sstart, sub.end as send,
       gene.start as gstart, gene.end as gend
from mite
    join genome on genome.id = mite.genome_id
    join gene on gene.genome_id = mite.genome_id and
         gene.name = mite.correlated_gene
    join gene_sub_feature as sub on sub.gene_name = gene.name and
         sub.genome_id = gene.genome_id limit 10;

+-----------+------+----------+----------+-----------------+----------+----------+----------+----------+
| name      | code | start    | end      | type            | sstart   | send     | gstart   | gend     |
+-----------+------+----------+----------+-----------------+----------+----------+----------+----------+
| AT1M10000 | sf   | 13999039 | 13999549 | exon            | 13998637 | 13999051 | 13998637 | 13999051 |
| AT1M10000 | sf   | 13999039 | 13999549 | five_prime_UTR  | 13999041 | 13999051 | 13998637 | 13999051 |
| AT1M10000 | sf   | 13999039 | 13999549 | three_prime_UTR | 13998637 | 13998843 | 13998637 | 13999051 |
| AT1M10010 | col  | 13944996 | 13945397 | exon            | 13944917 | 13945581 | 13944917 | 13945581 |
| AT1M10010 | col  | 13944996 | 13945397 | five_prime_UTR  | 13945571 | 13945581 | 13944917 | 13945581 |
| AT1M10010 | col  | 13944996 | 13945397 | three_prime_UTR | 13944917 | 13945112 | 13944917 | 13945581 |
| AT1M10040 | mt   | 14010954 | 14011023 | exon            | 14010700 | 14011105 | 14010700 | 14011105 |
| AT1M10040 | mt   | 14010954 | 14011023 | five_prime_UTR  | 14011101 | 14011105 | 14010700 | 14011105 |
| AT1M10040 | mt   | 14010954 | 14011023 | three_prime_UTR | 14010700 | 14010888 | 14010700 | 14011105 |
| AT1M10050 | mt   | 14010985 | 14011023 | exon            | 14010700 | 14011105 | 14010700 | 14011105 |
+-----------+------+----------+----------+-----------------+----------+----------+----------+----------+
```
