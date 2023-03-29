from .create_task import CreateTask
from .solve_task import SolveTask
from .free_mode import FreeMode
from ...bases.menu_tool import MenuTool

menu_tools: tuple[type[MenuTool], ...] = (
    CreateTask,
    SolveTask,
    FreeMode,
)
