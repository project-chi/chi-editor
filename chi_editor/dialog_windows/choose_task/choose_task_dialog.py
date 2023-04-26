from typing import TYPE_CHECKING
from random import choice

from PyQt6.QtWidgets import QDialog, QTreeView, QSizePolicy, QVBoxLayout, QAbstractItemView
from PyQt6.QtGui import QStandardItemModel, QStandardItem
from PyQt6.QtCore import Qt, QModelIndex

from chi_editor.api.task import Task, Kind

from chi_editor.editor_mode import EditorMode

if TYPE_CHECKING:
    from chi_editor.editor import Editor


class ChooseTaskDialog(QDialog):
    # Main window
    editor: "Editor"

    # View that holds all the tasks links
    view: QTreeView

    # Layout that holds view to make it expandable
    layout: QVBoxLayout

    # Model that links to all the tasks
    model: QStandardItemModel

    # Mapping from kinds to their entries in model
    kind_items: dict[Kind, QStandardItem]

    def __init__(self, *args, editor: "Editor", **kwargs) -> None:
        super().__init__(*args, **kwargs)
        self.setWindowTitle("Choose a task")

        self.editor = editor

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

    def handleAcceptClick(self):
        self.handleDoubleClick(self.view.currentIndex())

    def handleDoubleClick(self, index: QModelIndex) -> None:
        task_item = self.model.itemFromIndex(index)
        task = task_item.data(Qt.ItemDataRole.UserRole)
        if isinstance(task, Task):
            self.chooseTask(task)
            self.editor.setFormulationOfTask()

    def chooseTask(self, task: Task) -> None:
        self.editor.setMode(EditorMode.SOLVE_MODE)
        self.editor.setTask(task=task)
        self.close()
