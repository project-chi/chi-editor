from typing import TYPE_CHECKING
from random import choice

from PyQt6.QtWidgets import QTreeView, QSizePolicy, QVBoxLayout, QHBoxLayout, QPushButton, QAbstractItemView
from PyQt6.QtGui import QStandardItemModel, QStandardItem
from PyQt6.QtCore import Qt, QModelIndex

from chi_editor.api.server import Server, default_url
from chi_editor.api.task import Task, Kind

from chi_editor.editor_mode import EditorMode

from chi_editor.dialog_windows.choose_task.choose_task_dialog import ChooseTaskDialog

if TYPE_CHECKING:
    from chi_editor.editor import Editor


class RemoteTaskDialog(ChooseTaskDialog):
    # Main window
    editor: "Editor"

    # View that holds all the tasks links
    view: QTreeView

    # View's layout to hold buttons on top of tasks
    view_layout: QHBoxLayout

    # Buttons to manipulate items
    accept_button: QPushButton
    load_tasks_button: QPushButton
    random_task_button: QPushButton

    # Layout that holds view to make it expandable
    layout: QVBoxLayout

    # Model that links to all the tasks
    model: QStandardItemModel

    # Mapping from kinds to their entries in model
    kind_items: dict[Kind, QStandardItem]

    # Tasks database server
    server: Server

    def __init__(self, *args, editor: "Editor", **kwargs) -> None:
        super().__init__(*args, editor=editor, **kwargs)

        self.server = Server(default_url)
        
        # View layout
        self.view_layout = QHBoxLayout(self.view)
        self.view_layout.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignBottom)
        self.view_layout.setContentsMargins(0, 0, 2, 2)

        # Buttons
        self.setButtons()

    def setButtons(self) -> None:
        self.load_tasks_button = QPushButton("Load tasks")
        self.load_tasks_button.setFixedSize(self.load_tasks_button.sizeHint())
        self.load_tasks_button.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed)
        self.load_tasks_button.clicked.connect(self.loadTasks)

        self.accept_button = QPushButton("Choose task")
        self.accept_button.setFixedSize(self.accept_button.sizeHint())  # sizeHint() is minimal size to fit the text
        self.accept_button.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed)
        self.accept_button.clicked.connect(self.handleAcceptClick)

        self.random_task_button = QPushButton("Get random task")
        self.random_task_button.setFixedSize(self.random_task_button.sizeHint())
        self.random_task_button.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed)
        self.random_task_button.clicked.connect(self.handleRandomTaskClick)

        self.view_layout.addWidget(self.load_tasks_button)
        self.view_layout.addWidget(self.accept_button)
        self.view_layout.addWidget(self.random_task_button)

    def loadTasks(self) -> None:
        self._clearTasksList()

        tasks = self.server.get_tasks_raw()

        for task in tasks:
            self.addTask(task)
