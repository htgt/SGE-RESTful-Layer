from typing import List


class BaseClass:
    def get_fields(self) -> List[str]:
        keys = []
        for key in vars(self).keys():
            keys.append(str(key))
        return list(set(keys))

    def _asdict(self) -> dict:
        return vars(self)


class BaseConnection(BaseClass):
    pass
    
