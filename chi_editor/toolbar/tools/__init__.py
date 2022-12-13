from .arrow import Arrow
from .eraser import Eraser
from .nitrogen import Nitrogen
from .oxygen import Oxygen
from .structure import Structure
from .text import Text
from .bond import Bond
from .atom import Atom
from ...bases.tool import Tool


tools: tuple[type[Tool], ...] = (
    Arrow,
    Text,
    Bond,
    Atom,
    Nitrogen,
    Oxygen,
    Structure,
    Eraser,
)
