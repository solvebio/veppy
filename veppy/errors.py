class VeppyFeatureException(Exception):
    pass


class StopEffectPrediction(Exception):
    pass


class VeppyFileException(Exception):
    pass


class FeatureFileException(VeppyFileException):
    pass


class FastaFileException(VeppyFileException):
    pass
