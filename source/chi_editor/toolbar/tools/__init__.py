from .arrow import Arrow
from .text import Text
from ...abc.tool import Tool


tools: tuple[type[Tool], ...] = (
    Arrow,
    Text,
)
