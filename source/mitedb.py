# Copyright (C) 2015 Matthew R. Spinelli
# Department of Biology, Miami University

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

from schema import *
from sqlalchemy import select
import re
from itertools import chain


GFF, CALC = 0, 1


def parse_mite_gene_nearby(file_path, gene_chromos):
    import csv
    data = dict()
    with open(file_path) as in_handle:
        csv_reader = csv.reader(in_handle, delimiter='\t')
        for gene, gstart, gend, motif, mstart, mend in csv_reader:
            if not gene.startswith('novel'):
                chromo = gene[2]
            else:
                chromo = gene_chromos[gene]
            data[gene, int(gstart)-1, int(gend)] = \
                motif, int(mstart)-1, int(mend)
            data[chromo, int(mstart)-1, int(mend)] = \
                gene, int(gstart)-1, int(gend)
    return data


def lookup_chromosome(gene_file):
    gene_chromos = dict()
    for feature in parse_gff(gene_file, limit_info=dict(gff_type=['mRNA'])):
        pattern = '(novel_gene_[0-9]+)'
        gene_match = re.search(pattern, feature.id)
        if gene_match:
            feature.id = gene_match.group(1)
            gene_chromos[feature.id] = feature.chromosome
    return gene_chromos


def parse_gff(file_path, limit_info=None):
    from BCBio import GFF
    with open(file_path) as in_handle:
        for rec in GFF.parse(in_handle, limit_info=limit_info):
            for feature in rec.features:
                feature.chromosome = rec.id.strip('Chr')
                yield feature


def parse_genes(file_path, has_mite_in_at_least_one_accession, limit_info=None):
    for feature in parse_gff(file_path, limit_info):
        if not feature.sub_features:
            continue  # We don't want to import unannotated records right?
                      # Also, I'm getting the same mRNA entry multiple times
                      # with the second having no annotation: sf, AT1G02470
        pattern = '([A-Z]{2,}.+?[A-Z]+[0-9]+)'
        gene_match = re.search(pattern, feature.id)
        if not gene_match:
            pattern = '(novel_gene_[0-9]+)'
            gene_match = re.search(pattern, feature.id)
        if not gene_match:
            message = 'Cannot parse id from {0}'
            raise ImportError(message.format(file_path))
        if gene_match:
            feature.id = gene_match.group(1)
            key = (feature.id, int(feature.location.start),
                   int(feature.location.end))
            if feature.id not in has_mite_in_at_least_one_accession:
                continue
        annotate_sub_feature_counts(feature)
        yield feature


def parse_mites(file_path, nearby, limit_info=None):
    i = 0
    for feature in parse_gff(file_path, limit_info):
        i += 1
        target = feature.qualifiers['Target']
        if not target:
            message = 'Cannot target from {0}'
            raise ImportError(message.format(file_path))
        match = re.search(r'\"Motif:(.+)\"', target[0])
        if not match:
            message = 'Cannot motif from {0}'
            raise ImportError(message.format(file_path))
        motif = match.group(1)
        feature.motif = motif
        chromosome = feature.chromosome
        key = chromosome, int(feature.location.start), int(feature.location.end)
        feature.correlated_gene = nearby[key][0] if key in nearby else None
        feature.id = 'AT{0}M{1}'.format(chromosome, str(i * 10))
        yield feature


def get_mite_gene_location(mite, gene_sub_features):
    start_feature, end_feature = None, None
    is_reverse_strand = gene_sub_features[0].strand == '-'
    gene_sub_features.sort(key=lambda x: x.start)
    introns = [[gene_sub_features[0].end, None]]
    start = mite.location.start
    end = mite.location.end
    for sub in gene_sub_features:
        if sub.start <= start <= sub.end:
            if not start_feature:
                start_feature = sub.name
            else:
                end_feature = sub.name
        if sub.start <= end <= sub.end:
            if not start_feature:
                start_feature = sub.name
            else:
                end_feature = sub.name
        # begin intron interpolation
        intron_start = introns[-1][0]
        if intron_start < sub.start:
            introns[-1][1] = sub.start - 1
            introns.append([sub.end, None])
        else:
            introns[-1] = [sub.end, None]
    introns = introns[:-1] # remove last (it's not needed)
    if (start_feature and end_feature) or not introns:
        return start_feature, end_feature
    # Must be in an intron somewhere, so see if mite starts or ends in introns
    intron_count = 0
    if is_reverse_strand:
        introns = introns[::-1]
    for intron in introns:
        intron_count += 1
        if intron[0] <= start <= intron[1]:
            if not start_feature:
                start_feature = 'I' + str(intron_count)
            else:
                end_feature = 'I' + str(intron_count)
        if intron[0] <= end <= intron[1]:
            if not start_feature:
                start_feature = 'I' + str(intron_count)
            else:
                end_feature = 'I' + str(intron_count)
    return start_feature, end_feature


def annotate_sub_feature_counts(feature):
    """
    Calculates the number of exons (including some 5' and 3' UTRs as exons if
    they are spliced. Assumption: One 5' and 3' UTR sub feature each are
    adjacent to an exon and thus two of them total are not counted.

    Also assigns number to utrs, exons, and introns in the order they appear in
    the sequence.
    """
    exons, introns, utr5s, utr3s = 0, 0, 0, 0
    reverse = feature.strand == -1
    feature.sub_features.sort(
        key= lambda x: x.location.start, reverse=reverse)
    for sub in feature.sub_features:
        if sub.type == 'five_prime_UTR':
            utr5s += 1
            sub.name = 'F' + str(utr5s)
        elif sub.type == 'three_prime_UTR':
            utr3s += 1
            sub.name = 'T' + str(utr3s)
        elif sub.type == 'exon':
            exons += 1
            sub.name = 'E' + str(exons)
        else:
            sub.name = None
    feature.exon_count = exons


def annotate_starting_ending_features(conn, mites, code):
    for mite in mites:
        start, end = None, None
        if mite.correlated_gene:
            gene_sub_features = get_gene(conn, mite, code)
            if gene_sub_features:
                start, end = get_mite_gene_location(mite, gene_sub_features)
        mite.starting_feature = start
        mite.ending_feature = end
        yield mite


def suppress_warnings():
    from warnings import filterwarnings
    import pymysql
    filterwarnings('ignore', category=pymysql.Warning)


def connect_to_db(sql_url, database, new=False, echo=False):
    from sqlalchemy import create_engine
    suppress_warnings()
    engine = create_engine(sql_url, echo=echo)
    if new:
        if 'sqlite' not in sql_url:
            engine.execute('CREATE DATABASE IF NOT EXISTS {0}'.format(database))
            engine.execute('USE {0}'.format(database))
        metadata.drop_all(engine)
        metadata.create_all(engine)
        return engine
    if 'sqlite' not in sql_url:
        engine.execute('USE {0}'.format(database))
    return engine


def get_gene(conn, mite, code):
    from sqlalchemy import text
    genome_id = lookup_genome_id_from_code(conn, code)
    stmt = '''
        SELECT sub.name, sub.start, sub.end, gene.strand
        FROM gene_sub_feature AS sub
        JOIN gene ON gene.name = sub.gene_name AND
                     gene.genome_id = sub.genome_id
        WHERE sub.gene_name = :mite AND sub.genome_id = :id
        '''
    return conn.execute(text(stmt), mite=mite.correlated_gene,
                        id=genome_id).fetchall()


def save_genome(conn, code, accession, taxonomy, comments = None):
    stmt = db_genome.insert().prefix_with('IGNORE').values({
        'code': code,
        'accession': accession,
        'taxonomy': taxonomy,
        'comments': comments})
    conn.execute(stmt)


def lookup_genome_id_from_code(conn, code):
    stmt = select([db_genome.c.id]).where(db_genome.c.code == code)
    genome_ids = conn.execute(stmt).fetchone()
    return genome_ids[0] if genome_ids else None


def save_genes(conn, genes, code, accession, taxonomy, comments=None):
    genome_id = lookup_genome_id_from_code(conn, code)
    if not genome_id:
        save_genome(conn, code, accession, taxonomy, comments)
        genome_id = lookup_genome_id_from_code(conn, code)
    gene_stmts = []
    gene_sub_feature_stmts = []
    for feature in genes:
        strand = feature.strand
        gene_stmts.append({
            'name': feature.id,
            'start': int(feature.location.start),
            'end': int(feature.location.end),
            'chromosome': feature.chromosome,
            'strand': '+' if strand == 1 else '-' if strand == -1 else '.',
            'exon_count': feature.exon_count,
            'genome_id': genome_id,
            'accession_id': code})
        for sub_feature in feature.sub_features:
            gene_sub_feature_stmts.append({
                'name': sub_feature.name,
                'gene_name': feature.id,
                'type': sub_feature.type,
                'start': int(sub_feature.location.start),
                'end': int(sub_feature.location.end),
                'genome_id': genome_id})
    if gene_stmts:
        conn.execute(db_gene.insert().prefix_with('IGNORE'), gene_stmts)
    if gene_sub_feature_stmts:
        conn.execute(db_gene_sub_feature.insert().prefix_with('IGNORE'),
                     gene_sub_feature_stmts)


def save_mites(conn, mites, code, accession, taxonomy, comments=None):
    genome_id = lookup_genome_id_from_code(conn, code)
    if not genome_id:
        save_genome(conn, code, accession, taxonomy, comments)
        genome_id = lookup_genome_id_from_code(conn, code)
    mite_stmts = []
    for feature in mites:
        strand = feature.strand
        parts = feature.motif.split('|')
        mite_stmts.append({
            'name': feature.id,
            'start': int(feature.location.start),
            'end': int(feature.location.end),
            'chromosome': feature.chromosome,
            'strand': '+' if strand == 1 else '-' if strand == -1 else '.',
            'correlated_gene': feature.correlated_gene,
            'starting_feature': feature.starting_feature,
            'ending_feature': feature.ending_feature,
            'superfamily': parts[-1],
            'family': parts[-2],
            'tir': parts[3],
            'tsd': parts[4],
            'sequence': None,
            'len': None,
            'genome_id': genome_id})
    if mite_stmts:
        conn.execute(db_mite.insert().prefix_with('IGNORE'), mite_stmts)


def populate_db_script():
    import sys
    import os

    # create database
    db_name = 'MITEdb'
    sql_url = 'mysql+pymysql://root:BioInfoLab@localhost'
    message = 'Creating {0} db (overriding if exists)... '.format('mysql')
    sys.stdout.write(message)
    sys.stdout.flush()
    engine = connect_to_db(sql_url, db_name, True, echo=False)
    conn = engine.connect()
    sys.stdout.write('Finished\n')
    sys.stdout.flush()

    # parse MITEs and add to the database from the calculated file
    path = os.path.join(os.getcwd(), 'mites')
    files = os.listdir(path)

    # build dictionary of accession numbers from the A.thaliana genome codes
    # Todo these aren't right, it's just a placeholder for later
    accession = dict(bur=-1, can=-2, col=-3, ct=-4, edi=-5, hi=-6, kn=-7,
                     ler=-8, mt=-9, no=-10, oy=-11, po=-12, rsch=-13, sf=-14,
                     tsu=-15, wil=-16, ws=-17, wu=-18, zu=-19)

    # build dictionary of file types
    genome_files = dict()
    for file in files:
        parts = file.split('.')
        code = parts[0]
        extension = parts[-1]
        if '~' in extension:
            continue
        if extension == 'bed':
            continue
        if code not in genome_files:
            genome_files[code] = [None, None]
        if extension == 'gff':
            i = GFF
        elif extension == 'CalculationOfMITEsAndGenes':
            i = CALC
        genome_files[code][i] = os.path.join(path, file)

    # now write to database and calculate
    nearby = {}  # dict[mite]=gene
    path = os.path.join(os.getcwd(), 'genes')
    gene_files = dict((x.split('.')[-2], x) for x in os.listdir(path))
    for code, files in genome_files.items():
        # build filter of nearby mites/genes
        # uncomment to speed up for testing
        # if code not in ['sf', 'col']: continue
        file_name = files[CALC].split('/')[-1]
        message = 'Reading MITE/Gene correlations from {0}... '
        sys.stdout.write(message.format(file_name))
        sys.stdout.flush()
        gene_file = os.path.join(path, gene_files[code])
        gene_chromos = lookup_chromosome(gene_file)
        data = parse_mite_gene_nearby(files[CALC], gene_chromos)
        nearby[code] = data
        print('Finished')

    # determine which genes have mites in at least one accession
    has_mite_in_at_least_one_accession = set(chain.from_iterable(
        (key[0] for key in data.keys()) for data in nearby.values()))

    # parse gene files and add to the database
    for file in gene_files.values():
        file_name = file.split('/')[-1]
        code = file.split('.')[-2]
        # uncomment to speed up for testing
        # if code not in ['sf', 'col']: continue
        limit_info = dict(gff_type=['mRNA', 'exon', 'five_prime_UTR',
                                    'three_prime_UTR'])
        sys.stdout.write('Importing Genes from {0}... '.format(file_name))
        sys.stdout.flush()
        genes = list(parse_genes(os.path.join(path, file),
                                 has_mite_in_at_least_one_accession,
                                 limit_info))
        save_genes(conn, genes, code, accession[code], taxonomy=3702)
        print('Finished')

    # parse mite files and add to the database
    for code, files in genome_files.items():
        # uncomment to speed up for testing
        # if code not in ['sf', 'col']: continue
        file_name = files[GFF].split('/')[-1]
        sys.stdout.write('Importing MITEs from {0}... '.format(file_name))
        sys.stdout.flush()
        mites = parse_mites(files[GFF], nearby[code])
        mites = annotate_starting_ending_features(conn, mites, code)
        save_mites(conn, mites, code, accession[code], taxonomy=3702)
        print('Finished')

if __name__ == '__main__':
    populate_db_script()
