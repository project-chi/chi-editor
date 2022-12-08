from .arrow import Arrow
from .text import Text
from ...bases.tool import Tool


tools: tuple[type[Tool], ...] = (
    Arrow,
    Text,
)
