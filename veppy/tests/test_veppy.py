import os
import sys
import logging
from itertools import islice
from unittest import TestCase

from veppy.veppy import calculate_consequences
from veppy.utils.readers import VcfReader

logger = logging.getLogger('test')


IGNORED_EFFECTS = {
    'feature_truncation',
    'feature_elongation',
    'coding_sequence_variant'
}


class VepVcfTest(TestCase):
    TEST_BUILD = 'GRCh37'
    GENE_MODEL = 'gencode'

    def get_absolute_filepath(self, filepath):
        return os.path.join(
            os.path.dirname(__file__),
            'data',
            os.path.basename(filepath)
        )

    def get_id_without_version(self, tid):
        return tid.split('.', 1)[0]

    def assertTest(self, test):
        (transcript_id, variant, required_effects, disallowed_effects) = test
        disallow_all = len(disallowed_effects) and disallowed_effects[0] == '*'
        results = calculate_consequences(
            self.TEST_BUILD,
            *variant,
            gene_model=self.GENE_MODEL
        ).results

        transcript_results = dict(
            [(self.get_id_without_version(r['transcript'].id),
              r['consequences'])
             for r in results
             ]
        )

        def _eff(csq):
            return csq.so_term

        def _eff_set(csqs):
            return set(map(_eff, csqs))

        try:
            # HACK!
            required_effects = set(required_effects) - IGNORED_EFFECTS

            # ensure there actually **are** required effects
            if required_effects:
                self.assertIn(transcript_id, transcript_results)

                for effect in required_effects:
                    self.assertIn(effect, _eff_set(transcript_results[transcript_id]))      # noqa

            if not disallow_all:
                for not_effect in disallowed_effects:
                    self.assertNotIn(not_effect, _eff_set(transcript_results[transcript_id]))      # noqa

            else:
                pruned_results = set([
                    r for r in _eff_set(transcript_results[transcript_id])
                    if r not in IGNORED_EFFECTS
                ])

                self.assertEqual(
                    set(required_effects), set(pruned_results)
                )

        except AssertionError, e:
            print '-------------------------------------'
            print '-------------------------------------'
            print '   VARIANT: ', str(variant)
            print '  REQUIRED: ', list(required_effects)
            print 'DISALLOWED: ', disallowed_effects
            print '   RESULTS: ', map(
                str,
                transcript_results.get(transcript_id, [])
            )
            print '-------------------------------------'
            print '-------------------------------------'
            raise AssertionError(e)

    def _test(self, filepath):
        def vcf_row_to_test(vcf_row):
            return (
                vcf_row['info']['TRANSCRIPT'][0],
                (
                    vcf_row['chromosome'],
                    vcf_row['start'],
                    vcf_row['reference_allele'],
                    vcf_row['alternate_alleles'][0]
                ),
                [x for x in vcf_row['info'].get('REQUIRED', []) if x],
                [x for x in vcf_row['info'].get('DISALLOWED', []) if x]
            )

        with VcfReader(self.get_absolute_filepath(filepath)) as reader:
            for row in islice(reader, sys.maxint):
                self.assertTest(vcf_row_to_test(row))

    def test_snps(self):
        self._test('snps.vcf')

    def test_insertions(self):
        self._test('insertions.vcf')

    def test_splice_regions(self):
        self._test('splice_regions.vcf')

    def test_deletions(self):
        self._test('deletions.vcf')

    def test_substitutions(self):
        self._test('substitutions.vcf')

    def test_difficults(self):
        self._test('difficults.vcf')

    # def test_difficults_chrX(self):
        # FIXME: Local does not have chrX
        # self._test('difficults.chrX.vcf')


class MostSevereEffectTestCase(TestCase):
    TEST_BUILD = 'GRCh37'
    GENE_MODEL = 'gencode'

    def test_single_gene(self):
        variant = ('1', 8025384, 'A', 'T')

        # most_severe = False (base). return results for multiple
        #  transcripts
        r0 = \
            calculate_consequences(
                self.TEST_BUILD,
                *variant
            ).results

        # most_severe = True. return most severe transcript
        r1 = calculate_consequences(
            self.TEST_BUILD,
            *variant,
            most_severe=True,
            default_only=True
        ).results

        self.assertGreater(len(r0), 1)
        self.assertEqual(len(r1), 1)

    def test_exclude_neighbors(self):
        # this variant has a neighboring gene, with exclude it should only
        # return its exact transcript (SF3B4) and not MTMR11 (the neighbor)
        variant = ('1', 149895762, 'G', 'T')

        r0 = \
            calculate_consequences(
                self.TEST_BUILD,
                *variant,
                exclude_neighbors=True
            ).results

        self.assertEqual(len(r0), 1)
        self.assertEqual(r0[0]['transcript'].gene.info['gene_name'],
                         'SF3B4')

        r1 = \
            calculate_consequences(
                self.TEST_BUILD,
                *variant
            ).results

        self.assertGreater(len(r1), 1)

    # def test_multi_gene(self):
    #     variant = ('3', 37089131, 'A', 'C')

    #     # most_severe = False (base). return results for multiple
    #     #  transcripts
    #     r0 = \
    #         calculate_consequences(self.TEST_BUILD, *variant).results

    #     # most_severe = True. return most severe transcript per gene...
    #     #  this variant overlaps transcripts from two separate genes
    #     r1 = calculate_consequences(
    #         self.TEST_BUILD, *variant, most_severe=True).results

    #     self.assertGreater(len(r0), 1)
    #     self.assertEqual(len(r1), 2)
