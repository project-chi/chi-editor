from .answer_field_tool import AnswerFieldTool
from ...bases.tool import Tool
from .arrow import Arrow
from .atoms.carbon import Carbon
from .atoms.nitrogen import Nitrogen
from .atoms.oxygen import Oxygen
from .bond import Bond
from .bonds.create_double_bond import CreateDoubleBond
from .bonds.create_single_bond import CreateSingleBond
from .bonds.create_triple_bond import CreateTripleBond
from .bonds.create_wedge_bond import CreateWedgeBond
from .drag import Drag
from .eraser import Eraser
from .smiles import Smiles
from .structure import Structure
from .text import Text

tools: tuple[type[Tool], ...] = (
    Arrow,
    Text,
    CreateSingleBond,
    CreateDoubleBond,
    CreateTripleBond,
    Carbon,
    Nitrogen,
    Oxygen,
    Structure,
    Smiles,
    Eraser,
    Drag,
    AnswerFieldTool,
)
