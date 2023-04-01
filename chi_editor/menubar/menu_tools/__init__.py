from ...bases.menu_tool import MenuTool
from .choose_task import ChooseTask
from .create_task import CreateTask
from .free_mode import FreeMode
from .solve_task import SolveTask

mode_tools: tuple[type[MenuTool], ...] = (
    CreateTask,
    # SolveTask,
    FreeMode,
)

task_tools: tuple[type[MenuTool], ...] = (
    ChooseTask,
)
