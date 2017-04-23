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

from mitedb import *
from Bio.SeqFeature import SeqFeature, FeatureLocation, ExactPosition
from unittest import TestCase
from collections import namedtuple


SubFeature = namedtuple('SubFeature', field_names='strand, start, end, name')


class TestMitedb(TestCase):
    def test_get_mite_gene_location_exon_reverse(self):
        # Setup
        mite = SeqFeature(FeatureLocation(ExactPosition(14301135-1),
                                          ExactPosition(14301495),
                                          strand=-1),
                          type='mRNA', id='AT1G38630')

        # note that this is real data (can accession code)
        sub_features = [
            SubFeature(strand='-', start=14298853-1, end=14299101, name='T1'),
            SubFeature(strand='-', start=14298853-1, end=14299175, name='E7'),
            SubFeature(strand='-', start=14299460-1, end=14299528, name='E6'),
            SubFeature(strand='-', start=14301089-1, end=14301157, name='E5'),
            SubFeature(strand='-', start=14301443-1, end=14301511, name='E4'),
            SubFeature(strand='-', start=14301621-1, end=14301689, name='E3'),
            SubFeature(strand='-', start=14302679-1, end=14302747, name='E2'),
            SubFeature(strand='-', start=14302843-1, end=14302939, name='E1'),
            SubFeature(strand='-', start=14302892-1, end=14302939, name='F1'),]

        expected_start, expected_end = 'E5', 'E4'

        # Exercise
        start, end = get_mite_gene_location(mite, sub_features)

        # Verify
        self.assertEqual(start, expected_start)
        self.assertEqual(end, expected_end)

    def test_get_mite_gene_location_exon_forward(self):
        # Setup
        mite = SeqFeature(FeatureLocation(ExactPosition(14301135-1),
                                          ExactPosition(14301495),
                                          strand=1),
                          type='mRNA', id='AT1G38630')

        sub_features = [
            SubFeature(strand='+', start=14298853-1, end=14299101, name='F1'),
            SubFeature(strand='+', start=14298853-1, end=14299175, name='E1'),
            SubFeature(strand='+', start=14299460-1, end=14299528, name='E2'),
            SubFeature(strand='+', start=14301089-1, end=14301157, name='E3'),
            SubFeature(strand='+', start=14301443-1, end=14301511, name='E4'),
            SubFeature(strand='+', start=14301621-1, end=14301689, name='E5'),
            SubFeature(strand='+', start=14302679-1, end=14302747, name='E6'),
            SubFeature(strand='+', start=14302843-1, end=14302939, name='E7'),
            SubFeature(strand='+', start=14302892-1, end=14302939, name='T1'),]

        expected_start, expected_end = 'E3', 'E4'

        # Exercise
        start, end = get_mite_gene_location(mite, sub_features)

        # Verify
        self.assertEqual(start, expected_start)
        self.assertEqual(end, expected_end)

    def test_get_mite_gene_location_intron_reverse(self):
        # Setup
        mite = SeqFeature(FeatureLocation(ExactPosition(9762301-1),
                                          ExactPosition(9762350),
                                          strand=-1),
                          type='mRNA', id='AT1G28230')

        sub_features = [
            SubFeature(strand='-', start=9761599-1, end=9761802, name='T1'),
            SubFeature(strand='-', start=9761599-1, end=9762165, name='E2'),
            # SubFeature(strand='-', start=9762166-1, end=9763449, name='I1'), this is what should be calculated
            SubFeature(strand='-', start=9763450-1, end=9764167, name='E1'),
            SubFeature(strand='-', start=9764158-1, end=9764167, name='F1'),]

        expected_start, expected_end = 'I1', 'I1'

        # Exercise
        start, end = get_mite_gene_location(mite, sub_features)

        # Verify
        self.assertEqual(start, expected_start)
        self.assertEqual(end, expected_end)

    def test_get_mite_gene_location_intron_forward(self):
        # Setup
        mite = SeqFeature(FeatureLocation(ExactPosition(9762301-1),
                                          ExactPosition(9762350),
                                          strand=1),
                          type='mRNA', id='AT1G28230')

        sub_features = [
            SubFeature(strand='+', start=9761599-1, end=9761802, name='T1'),
            SubFeature(strand='+', start=9761599-1, end=9762165, name='E2'),
            # SubFeature(strand='+', start=9762166-1, end=9763449, name='I1'), this is what should be calculated
            SubFeature(strand='+', start=9763450-1, end=9764167, name='E1'),
            SubFeature(strand='+', start=9764158-1, end=9764167, name='F1'),]

        expected_start, expected_end = 'I1', 'I1'

        # Exercise
        start, end = get_mite_gene_location(mite, sub_features)

        # Verify
        self.assertEqual(start, expected_start)
        self.assertEqual(end, expected_end)

    def test_get_mite_gene_location_intron_reverse_lots_of_introns(self):
        # Setup
        mite = SeqFeature(FeatureLocation(ExactPosition(511777-1),
                                          ExactPosition(512242),
                                          strand=1),
                          type='mRNA', id='AT1G02470')

        sub_features = [
            SubFeature(strand='-', start=510853-1, end=511011, name='T1'),
            SubFeature(strand='-', start=510853-1, end=511086, name='E7'),
            SubFeature(strand='-', start=511170-1, end=511217, name='E6'),
            SubFeature(strand='-', start=511310-1, end=511358, name='E5'),
            SubFeature(strand='-', start=511474-1, end=511526, name='E4'),
            SubFeature(strand='-', start=511621-1, end=511716, name='E3'),
            SubFeature(strand='-', start=512243-1, end=512342, name='E2'),
            SubFeature(strand='-', start=512428-1, end=512707, name='E1'),
            SubFeature(strand='-', start=512670-1, end=512707, name='F1'),]

        expected_start, expected_end = 'E2', 'I2'

        # Exercise
        start, end = get_mite_gene_location(mite, sub_features)

        # Verify
        self.assertEqual(start, expected_start)
        self.assertEqual(end, expected_end)


    def test_annotate_sub_feature_reverse_correct_annotation_and_counts(self):
        # Setup
        sub_features = [SeqFeature(FeatureLocation(ExactPosition(9761599-1),
                                                   ExactPosition(9761802),
                                                   strand=-1),
                                   type='three_prime_UTR'),
                        SeqFeature(FeatureLocation(ExactPosition(9761599-1),
                                                   ExactPosition(9762165),
                                                   strand=-1),
                                   type='exon'),
                        SeqFeature(FeatureLocation(ExactPosition(9763450-1),
                                                   ExactPosition(9764167),
                                                   strand=-1),
                                   type='exon'),
                        SeqFeature(FeatureLocation(ExactPosition(9764158-1),
                                                   ExactPosition(9764167),
                                                   strand=-1),
                                   type='five_prime_UTR'),]

        feature = SeqFeature(FeatureLocation(ExactPosition(9762301-1),
                                             ExactPosition(9762350),
                                             strand=-1),
                             type='mRNA', id='AT1G28230',
                             sub_features=sub_features)

        expected_exon_counts = 2
        names = iter(['F1', 'E1', 'T1', 'E2'])  # possible sort diff if tie?

        # Exercise
        annotate_sub_feature_counts(feature)

        # Verify
        self.assertEqual(feature.exon_count, expected_exon_counts)
        for sub in feature.sub_features:
            self.assertEqual(sub.name, next(names))
