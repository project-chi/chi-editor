from PyQt6.QtWidgets import QDialog, QTreeView, QSizePolicy, QVBoxLayout, QHBoxLayout, QPushButton, QAbstractItemView
from PyQt6.QtGui import QStandardItemModel, QStandardItem
from PyQt6.QtCore import Qt, QModelIndex

from .tasks.task import Task
from .tasks import tasks_list


class ChooseTaskDialog(QDialog):
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

    def __init__(self, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs)
        self.setWindowTitle("Choose a task")

        # Model init
        self.model = QStandardItemModel()
        self.fillModel()

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
        self.accept_button = QPushButton("Choose task")
        self.accept_button.setFixedSize(self.accept_button.sizeHint())  # sizeHint() is minimal size to fit the text
        self.accept_button.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed)

        self.accept_button.clicked.connect(self.handleAcceptClick)

        self.view_layout.addWidget(self.accept_button)

    def fillModel(self) -> None:
        for tasks_desc in tasks_list:
            type_item = QStandardItem(tasks_desc[1].name)
            self.model.appendRow(type_item)
            for task in tasks_desc[0]:
                task_item = QStandardItem(task.title())
                task_item.setData(task, Qt.ItemDataRole.UserRole)  # UserRole means application specific role
                type_item.appendRow(task_item)

    def handleAcceptClick(self):
        current_index = self.view.currentIndex()
        task_item = self.model.itemFromIndex(current_index)
        task = task_item.data(Qt.ItemDataRole.UserRole)
        self.chooseTask(task)

    def handleDoubleClick(self, index: QModelIndex) -> None:
        task_item = self.model.itemFromIndex(index)
        task = task_item.data(Qt.ItemDataRole.UserRole)
        self.chooseTask(task)

    def chooseTask(self, task: Task) -> None:
        print(task.title())
        self.close()
