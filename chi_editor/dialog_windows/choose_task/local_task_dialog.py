from typing import TYPE_CHECKING, ClassVar
from pathlib import Path
from random import choice

from PyQt6.QtWidgets import QTreeView, QSizePolicy, QVBoxLayout, QHBoxLayout, QPushButton
from PyQt6.QtGui import QStandardItemModel, QStandardItem
from PyQt6.QtCore import Qt, QModelIndex

from chi_editor.constants import RESOURCES
from chi_editor.api.task import Task, Kind

from chi_editor.editor_mode import EditorMode

from chi_editor.dialog_windows.choose_task.choose_task_dialog import ChooseTaskDialog

if TYPE_CHECKING:
    from chi_editor.editor import Editor


class LocalTaskDialog(ChooseTaskDialog):
    # Default folder for local task files
    default_dir: ClassVar[Path] = RESOURCES / "local_tasks"

    def __init__(self, *args, editor: "Editor", **kwargs) -> None:
        super().__init__(*args, editor=editor, **kwargs)

    def loadTasks(self) -> None:
        self._clearTasksList()

        for json_file in self.default_dir.glob("*.json"):
            task = Task.parse_file(json_file)
            self.addTask(task)
