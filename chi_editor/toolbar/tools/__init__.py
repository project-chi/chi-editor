from ...bases.tool import Tool
from .arrow import Arrow
from .drag import Drag
from .eraser import Eraser
from .smiles import Smiles
from .structure import Structure
from .text import Text

from ...bases.toolbar_menu_widget import ToolbarMenuWidget
from .atoms_menu import AtomsMenu
from .atoms import atom_tools
from .bonds_menu import BondsMenu
from .bonds import bond_tools

tools: tuple[type[Tool], ...] = (
    Arrow,
    Text,
    Structure,
    Smiles,
    Eraser,
    Drag,
)

menus: list[tuple[tuple[type[Tool], ...], type[ToolbarMenuWidget]], ...] = [
    (atom_tools, AtomsMenu),
    (bond_tools, BondsMenu),
]
