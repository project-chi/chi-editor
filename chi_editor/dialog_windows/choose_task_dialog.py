from typing import TYPE_CHECKING
from random import choice

from PyQt6.QtWidgets import QDialog, QTreeView, QSizePolicy, QVBoxLayout, QHBoxLayout, QPushButton, QAbstractItemView
from PyQt6.QtGui import QStandardItemModel, QStandardItem
from PyQt6.QtCore import Qt, QModelIndex

from chi_editor.api.server import Server, default_url
from chi_editor.api.task import Task, Kind

from chi_editor.editor_mode import EditorMode

if TYPE_CHECKING:
    from chi_editor.editor import Editor


class ChooseTaskDialog(QDialog):
    # Main window
    editor: "Editor"

    # View that holds all the tasks links
    view: QTreeView

    # View's layout to hold buttons on top of tasks
    view_layout: QHBoxLayout

    # Buttons to manipulate items
    accept_button: QPushButton

    # Layout that holds view to make it expandable
    layout: QVBoxLayout

    # Model that links to all the tasks
    model: QStandardItemModel

    # Mapping from kinds to their entries in model
    kind_items: dict[Kind, QStandardItem]

    # Tasks database server
    server: Server

    def __init__(self, *args, editor: "Editor", **kwargs) -> None:
        super().__init__(*args, **kwargs)
        self.setWindowTitle("Choose a task")

        self.editor = editor
        self.server = Server(default_url)

        # Model init
        self.model = QStandardItemModel()
        self.kind_items = {}
        self.fillModelKinds()

        # View init
        self.view = QTreeView(self)
        self.view.setSizePolicy(QSizePolicy.Policy.MinimumExpanding, QSizePolicy.Policy.MinimumExpanding)
        self.view.setModel(self.model)
        self.view.doubleClicked.connect(self.handleDoubleClick)
        self.view.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.view.setHeaderHidden(True)

        # Main layout
        self.layout = QVBoxLayout(self)
        self.layout.setContentsMargins(0, 2, 0, 0)
        self.layout.addWidget(self.view)

        # View layout
        self.view_layout = QHBoxLayout(self.view)
        self.view_layout.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignBottom)
        self.view_layout.setContentsMargins(0, 0, 2, 2)

        # Buttons
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

    def fillModelKinds(self) -> None:
        for kind in Kind:
            type_item = QStandardItem(kind.name)
            type_item.setData(kind, Qt.ItemDataRole.UserRole)
            self.model.appendRow(type_item)
            self.kind_items.update({kind: type_item})

    def _clearTasksList(self) -> None:
        for r in range(0, self.model.rowCount()):  # run through top level categories and remove their contents (rows)
            kind_row = self.model.item(r)
            kind_row.removeRows(0, kind_row.rowCount())

    def loadTasks(self) -> None:
        self._clearTasksList()

        tasks = self.server.get_tasks_raw()

        for task in tasks:
            task_item = QStandardItem(task.name)
            task_item.setData(task, Qt.ItemDataRole.UserRole)

            kind_item = self.kind_items.get(task.kind)  # get item containing corresponding kind with dictionary
            kind_item.appendRow(task_item)

    def handleAcceptClick(self):
        self.handleDoubleClick(self.view.currentIndex())

    def handleDoubleClick(self, index: QModelIndex) -> None:
        task_item = self.model.itemFromIndex(index)
        task = task_item.data(Qt.ItemDataRole.UserRole)
        if isinstance(task, Task):
            self.chooseTask(task)
            self.editor.setFormulationOfTask()

    def handleRandomTaskClick(self) -> None:
        task_ids = self.server.get_tasks()
        task = self.server.get_task(choice(task_ids))
        self.chooseTask(task)
        self.editor.setFormulationOfTask()

    def chooseTask(self, task: Task) -> None:
        self.editor.setMode(EditorMode.SOLVE_MODE)
        self.editor.setTask(task=task)
        self.close()
