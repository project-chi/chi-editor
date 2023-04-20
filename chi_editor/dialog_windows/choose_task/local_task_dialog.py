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
    # Main window
    editor: "Editor"

    # View that holds all the tasks links
    view: QTreeView

    # View's layout to hold buttons on top of tasks
    view_layout: QHBoxLayout

    # Buttons to manipulate items
    accept_button: QPushButton
    random_task_button: QPushButton
    delete_button: QPushButton

    # Layout that holds view to make it expandable
    layout: QVBoxLayout

    # Model that links to all the tasks
    model: QStandardItemModel

    # Mapping from kinds to their entries in model
    kind_items: dict[Kind, QStandardItem]

    # Default folder for local task files
    default_dir: ClassVar[Path] = RESOURCES / "local_tasks"

    def __init__(self, *args, editor: "Editor", **kwargs) -> None:
        super().__init__(*args, editor=editor, **kwargs)

        # Fill model
        self._fillModel()

        # View layout
        self.view_layout = QHBoxLayout(self.view)
        self.view_layout.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignBottom)
        self.view_layout.setContentsMargins(0, 0, 2, 2)

        # Buttons
        self.setButtons()

    def _fillModel(self) -> None:
        for json_file in self.default_dir.glob("*.json"):
            task = Task.parse_file(json_file)
            task_item = QStandardItem(task.name)
            task_item.setData(task, Qt.ItemDataRole.UserRole)

            kind_item = self.kind_items.get(task.kind)  # get item containing corresponding kind with dictionary
            kind_item.appendRow(task_item)

    def setButtons(self) -> None:
        self.accept_button = QPushButton("Choose task")
        self.accept_button.setFixedSize(self.accept_button.sizeHint())  # sizeHint() is minimal size to fit the text
        self.accept_button.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed)
        self.accept_button.clicked.connect(self.handleAcceptClick)

        self.delete_button = QPushButton("Delete task")
        self.delete_button.setFixedSize(self.delete_button.sizeHint())  # sizeHint() is minimal size to fit the text
        self.delete_button.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed)
        self.delete_button.clicked.connect(self.handleDeleteClick)

        self.random_task_button = QPushButton("Get random task")
        self.random_task_button.setFixedSize(self.random_task_button.sizeHint())
        self.random_task_button.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed)
        self.random_task_button.clicked.connect(self.handleRandomTaskClick)

        self.view_layout.addWidget(self.accept_button)
        self.view_layout.addWidget(self.random_task_button)
        self.view_layout.addWidget(self.delete_button)

    def handleRandomTaskClick(self) -> None:
        pass
