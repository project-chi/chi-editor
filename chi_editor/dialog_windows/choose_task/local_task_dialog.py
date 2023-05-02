from typing import TYPE_CHECKING, ClassVar
from pathlib import Path

from PyQt6.QtWidgets import QTreeView, QSizePolicy, QDialog
from PyQt6.QtGui import QFileSystemModel

from chi_editor.constants import RESOURCES, ASSETS

from chi_editor.api.task import Task, Kind

from chi_editor.dialog_windows.choose_task.choose_task_dialog import ChooseTaskDialog

if TYPE_CHECKING:
    from chi_editor.editor import Editor


class LocalTaskDialog(ChooseTaskDialog):
    # Default folder for local task files
    default_dir: ClassVar[Path] = RESOURCES / "local_tasks"

    # Custom task directory
    task_dir: Path = default_dir

    # Local file system view/model
    dir_dialog: QDialog
    dir_view: QTreeView
    dir_model: QFileSystemModel

    def __init__(self, *args, editor: "Editor", **kwargs) -> None:
        super().__init__(*args, editor=editor, **kwargs)

        # Local file system
        self.dir_dialog = QDialog(self)
        self.dir_dialog.setWindowTitle("Choose directory")

        self.dir_view = QTreeView(self.dir_dialog)
        self.dir_view.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed)

        self.dir_dialog.setSizePolicy(QSizePolicy.Policy.MinimumExpanding, QSizePolicy.Policy.MinimumExpanding)

        self.dir_model = QFileSystemModel(self.dir_view)
        root_index = self.dir_model.setRootPath(str(self.default_dir.parent))
        self.dir_model.setReadOnly(True)

        self.dir_view.setModel(self.dir_model)
        self.dir_view.setRootIndex(root_index)

    def setSettingsActions(self) -> None:
        self.settings_menu.addAction("Change task directory", self.showDirChangeDialog)

    def showDirChangeDialog(self):
        self.dir_dialog.exec()

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
