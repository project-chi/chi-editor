from ...bases.tool import Tool
from .arrow import Arrow
from .drag import Drag
from .eraser import Eraser
from .smiles import Smiles
from .structure import Structure
from .text import Text

from .atoms_menu import AtomsMenu
from .atoms import atom_tools

tools: tuple[type[Tool], ...] = (
    Arrow,
    Text,
    Structure,
    Smiles,
    Eraser,
    Drag,
)

menus: list[tuple[tuple[type[Tool], ...], type[Tool]], ...] = [
    (atom_tools, AtomsMenu),
]
