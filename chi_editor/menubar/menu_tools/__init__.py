from .create_task import CreateTask
from .solve_task import SolveTask
from .free_mode import FreeMode
from ...bases.menu_tool import MenuTool

mode_tools: tuple[type[MenuTool], ...] = (
    CreateTask,
    SolveTask,
    FreeMode,
)

task_tools: tuple[type[MenuTool], ...] = (

)
