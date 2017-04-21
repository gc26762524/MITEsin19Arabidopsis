# Copyright (C) 2015 Matthew R. Spinelli
# Department of Biology, Miami University

from sqlalchemy import Table, MetaData, Column, Integer, String, ForeignKey

metadata = MetaData()

db_genome = Table(
    'genome', metadata,
    Column('id', Integer, autoincrement=True, primary_key=True),
    Column('accession', Integer, nullable=False),
    Column('taxonomy', Integer, nullable=False),
    Column('code', String(4), nullable=False),
    Column('comments', String(100), nullable=True),
)

db_mite = Table(
    'mite', metadata,
    Column('name', String(60), primary_key=True, nullable=False),
    Column('chromosome', String(20), nullable=False),
    Column('strand', String(1), nullable=False),
    Column('start', Integer, nullable=False),
    Column('end', Integer, nullable=False),
    Column('correlated_gene', String(60), nullable=True),
    Column('starting_feature', String(25), nullable=False),
    Column('ending_feature', String(25), nullable=False),
    Column('tir', Integer, nullable=False),
    Column('tsd', Integer, nullable=False),
    Column('family', String(20), nullable=False),
    Column('superfamily', String(20), nullable=False),
    Column('sequence', String(4000), nullable=False),
    Column('genome_id', Integer, ForeignKey(db_genome.c.id),
           primary_key=True, nullable=False)
)

db_gene = Table(
    'gene', metadata,
    Column('name', String(60), primary_key=True),
    Column('chromosome', String(20), nullable=False),
    Column('strand', String(1), nullable=False),
    Column('start', Integer, nullable=False),
    Column('end', Integer, nullable=False),
    Column('exon_count', Integer, nullable=False),
    Column('genome_id', Integer, ForeignKey(db_genome.c.id), primary_key=True,
           nullable=False)
)

db_gene_sub_feature = Table(
    'gene_sub_feature', metadata,
    Column('name', String(8), primary_key=True),
    Column('gene_name', String(60), ForeignKey(db_gene.c.name),
           primary_key= True),
    Column('type', String(20), nullable=False),
    Column('start', Integer, nullable=False),
    Column('end', Integer, nullable=False),
    Column('genome_id', Integer, ForeignKey(db_gene.c.genome_id),
           primary_key=True)
)