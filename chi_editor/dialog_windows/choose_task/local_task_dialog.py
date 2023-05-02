from typing import TYPE_CHECKING, ClassVar
from pathlib import Path

from PyQt6.QtWidgets import QTreeView, QSizePolicy, QDialog, QVBoxLayout, QHeaderView, QFileDialog
from PyQt6.QtGui import QFileSystemModel
from PyQt6.QtCore import QDir

from chi_editor.constants import RESOURCES

from chi_editor.api.task import Task

from chi_editor.dialog_windows.choose_task.choose_task_dialog import ChooseTaskDialog

if TYPE_CHECKING:
    from chi_editor.editor import Editor


class LocalTaskDialog(ChooseTaskDialog):
    # Default folder for local task files
    default_dir: ClassVar[Path] = RESOURCES / "local_tasks"

    # Layout that holds view to make it expandable
    main_layout: QVBoxLayout

    # Custom task directory
    task_dir: Path = default_dir

    # Change local directory
    dir_dialog: QDialog

    def __init__(self, *args, editor: "Editor", **kwargs) -> None:
        super().__init__(*args, editor=editor, **kwargs)

        # Local file system
        self.dir_dialog = QFileDialog(self)
        self.dir_dialog.setWindowTitle("Choose directory")

    def setSettingsActions(self) -> None:
        self.settings_menu.addAction("Change task directory", self.showDirChangeDialog)

    def showDirChangeDialog(self) -> None:
        self.dir_dialog.exec()

    def _addDirBottomButtons(self) -> None:
        pass

    def loadTasks(self) -> None:
        self._clearTasksList()

        for json_file in self.task_dir.glob("*.json"):
            task = Task.parse_file(json_file)
            self.addTask(task)

    def deleteTaskFromDatabase(self, task: Task) -> None:
        for file in self.task_dir.glob(task.name + ".json"):
            file.unlink()

    def updateWorkingDir(self, new_dir: Path) -> None:
        self.task_dir = new_dir
