from .create_task import CreateTask
from .solve_task import SolveTask
from .free_mode import FreeMode
from chi_editor.menubar.menu_tools.tasks.choose_remote_task import ChooseRemoteTask
from chi_editor.menubar.menu_tools.tasks.clear_tasks import ClearTasks
from ...bases.menu_tool import MenuTool

mode_tools: tuple[type[MenuTool], ...] = (
    CreateTask,
    # SolveTask,
    FreeMode,
)

task_tools: tuple[type[MenuTool], ...] = (
    ChooseRemoteTask,
    ClearTasks,
)
