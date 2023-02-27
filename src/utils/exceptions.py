class NoSecretKeyException(Exception):
    pass

class OligoDirectionInvalid(Exception):
    pass

### Warnings
class NoBenchlingEnvMatchWarning(UserWarning):
    pass