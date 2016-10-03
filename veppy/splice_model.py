from collections import namedtuple

FeatureRange = namedtuple('FeatureRange', ['start', 'stop'])


class BaseSpliceJunctionModel(object):
    @classmethod
    def get_acceptor_junction_model(klass, transcript, junction_model):
        return klass.AcceptorJunctionModel(transcript, junction_model)

    @classmethod
    def get_donor_junction_model(klass, transcript, junction_model):
        return klass.DonorJunctionModel(transcript, junction_model)


# SpliceJunctionModel is identical to VEP's
class SpliceJunctionModel(BaseSpliceJunctionModel):
    class AcceptorJunctionModel(object):
        def __init__(self, transcript, splice_junction_model):
            self.strand = transcript.strand
            self.model = splice_junction_model

            self.splice_acceptor = \
                FeatureRange(
                    self.model.splice_acceptor.start,
                    self.model.splice_acceptor.stop
                )
            if self.strand == '+':
                self.intron_splice_region = \
                    FeatureRange(
                        self.splice_acceptor.start - 6,
                        self.splice_acceptor.start - 1
                    )
                self.exon_splice_region = \
                    FeatureRange(
                        self.splice_acceptor.stop + 1,
                        self.splice_acceptor.stop + 3
                    )
            else:
                self.intron_splice_region = \
                    FeatureRange(
                        self.splice_acceptor.stop + 1,
                        self.splice_acceptor.stop + 6
                    )
                self.exon_splice_region = \
                    FeatureRange(
                        self.splice_acceptor.start - 3,
                        self.splice_acceptor.start - 1
                    )

    class DonorJunctionModel(object):
        def __init__(self, transcript, splice_junction_model):
            self.strand = transcript.strand
            self.model = splice_junction_model

            self.splice_donor = \
                FeatureRange(
                    self.model.splice_donor.start,
                    self.model.splice_donor.stop
                )

            if self.strand == '+':
                self.intron_splice_region = \
                    FeatureRange(
                        self.splice_donor.stop + 1,
                        self.splice_donor.stop + 6,
                    )
                self.exon_splice_region = \
                    FeatureRange(
                        self.splice_donor.start - 3,
                        self.splice_donor.start - 1
                    )
            else:
                self.intron_splice_region = \
                    FeatureRange(
                        self.splice_donor.start - 6,
                        self.splice_donor.start - 1,
                    )
                self.exon_splice_region = \
                    FeatureRange(
                        self.splice_donor.stop + 1,
                        self.splice_donor.stop + 3
                    )
