from typing import List
from dataclasses import dataclass

class BaseClass:
    def get_fields(self) -> List[str]:
        keys = []
        for key in vars(self).keys():
            keys.append(str(key))
        return list(set(keys))
    
@dataclass
class BaseBenchlingClass(BaseClass):
    targeton: str
    folder_id: str
    schema_id: str
    name: str
    strand: str
    grna: str