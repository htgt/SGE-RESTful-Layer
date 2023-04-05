class NoSecretKeyException(Exception):
    pass


class OligoDirectionInvalid(Exception):
    pass


class ConvertToJsonError(Exception):
    pass


class NoDotENVFile(Exception):
    pass

class MissingEnvVariables(Exception):
    pass

# Warnings


class NoBenchlingEnvMatchWarning(UserWarning):
    pass
