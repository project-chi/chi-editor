from typing import ClassVar
from pathlib import Path

from PyQt6.QtWidgets import QVBoxLayout, QFileDialog

from chi_editor.constants import RESOURCES

from chi_editor.api.task import Task

from chi_editor.dialog_windows.choose_task.choose_task_dialog import ChooseTaskDialog


class LocalTaskDialog(ChooseTaskDialog):
    # Default folder for local task files
    default_dir: ClassVar[Path] = RESOURCES / "local_tasks"

    # Layout that holds view to make it expandable
    main_layout: QVBoxLayout

    # Custom task directory
    task_dir: Path = default_dir

    def setSettingsActions(self) -> None:
        self.settings_menu.addAction("Change task directory", self.chooseLocalDirectory)

    def chooseLocalDirectory(self) -> None:
        new_local_dir_path = QFileDialog.getExistingDirectory(self, "Choose directory", str(self.task_dir.parent),
                                                              QFileDialog.Option.ShowDirsOnly)
        if new_local_dir_path == "":
            return
        self.task_dir = Path(new_local_dir_path)

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
