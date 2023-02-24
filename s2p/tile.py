from dataclasses import dataclass
from typing import List, Tuple


@dataclass(frozen=True)
class Tile:
    coordinates: Tuple[int, int, int, int]
    dir: str
    neighborhood_dirs: List[str]
    json: str
