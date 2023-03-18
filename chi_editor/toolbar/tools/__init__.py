from .arrow import Arrow
from .bonds.create_double_bond import CreateDoubleBond
from .bonds.create_single_bond import CreateSingleBond
from .bonds.create_triple_bond import CreateTripleBond
from .bonds.create_wedge_bond import CreateWedgeBond
from .eraser import Eraser
from .atoms.carbon import Carbon
from .atoms.nitrogen import Nitrogen
from .atoms.oxygen import Oxygen
from .structure import Structure
from .text import Text
from .bond import Bond
from .drag import Drag
from .smiles import Smiles
from ...bases.tool import Tool


tools: tuple[type[Tool], ...] = (
    Arrow,
    Text,
    CreateSingleBond,
    CreateDoubleBond,
    CreateTripleBond,
    CreateWedgeBond,
    Carbon,
    Nitrogen,
    Oxygen,
    Structure,
    Smiles,
    Eraser,
    Drag
)
