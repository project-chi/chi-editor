from typing import TYPE_CHECKING, cast
from random import randint

from PyQt6.QtWidgets import QDialog, QTreeView, QSizePolicy, QVBoxLayout, QAbstractItemView, QHBoxLayout, QPushButton, \
    QLineEdit, QMenu
from PyQt6.QtGui import QStandardItemModel, QStandardItem, QIcon
from PyQt6.QtCore import Qt, QModelIndex, QSortFilterProxyModel, QPoint

from chi_editor.api.task import Task, Kind
from chi_editor.editor_mode import EditorMode
from chi_editor.constants import ASSETS

if TYPE_CHECKING:
    from chi_editor.editor import Editor


class ChooseTaskDialog(QDialog):
    # Main window
    editor: "Editor"

    # View that holds all the tasks links
    view: QTreeView

    # Buttons to manipulate tasks
    accept_button: QPushButton
    random_task_button: QPushButton
    delete_button: QPushButton

    # Extra buttons
    load_tasks_button: QPushButton
    settings_button: QPushButton

    searchbar: QLineEdit

    # Layout that holds view to make it expandable
    layout: QVBoxLayout

    # Model that links to all the tasks
    model: QStandardItemModel
    proxy_model: QSortFilterProxyModel

    # Mapping from kinds to their entries in model
    kind_items: dict[Kind, QStandardItem]

    # Settings menu
    settings_menu: QMenu

    def __init__(self, *args, editor: "Editor", **kwargs) -> None:
        super().__init__(*args, **kwargs)
        self.setWindowTitle("Choose a task")

        self.editor = editor
        self.settings_menu = QMenu(self)

        # Model init
        self.model = QStandardItemModel()
        self.proxy_model = QSortFilterProxyModel()
        self.proxy_model.setSourceModel(self.model)
        self.proxy_model.setRecursiveFilteringEnabled(True)
        self.proxy_model.setFilterCaseSensitivity(Qt.CaseSensitivity.CaseInsensitive)
        self.proxy_model.setAutoAcceptChildRows(True)
        self.kind_items = {}

        # View init
        self.view = QTreeView(self)
        self.view.setModel(self.proxy_model)
        self.view.doubleClicked.connect(self.handleDoubleClick)
        self.view.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.view.setHeaderHidden(True)

        # Main layout
        self.layout = QVBoxLayout(self)
        self.layout.setContentsMargins(0, 2, 0, 0)

        self.setHeaderBar()
        self.layout.addWidget(self.view)
        self.setMainButtons()

        # Size configuration
        self.view.setMinimumSize(1, self.view.fontMetrics().height())
        self.view.setSizePolicy(QSizePolicy.Policy.MinimumExpanding, QSizePolicy.Policy.MinimumExpanding)

    def setMainButtons(self) -> None:
        # Buttons layout
        view_layout = QHBoxLayout(self)
        view_layout.setAlignment(Qt.AlignmentFlag.AlignRight)

        self.accept_button = QPushButton("Choose task")
        self.accept_button.clicked.connect(self.handleAcceptClick)

        self.random_task_button = QPushButton("Get random task")
        self.random_task_button.clicked.connect(self.handleRandomTaskClick)

        self.delete_button = QPushButton("Delete task")
        self.delete_button.clicked.connect(self.handleDeleteClick)

        view_layout.addWidget(self.accept_button)
        view_layout.addWidget(self.random_task_button)
        view_layout.addWidget(self.delete_button)

        self.layout.addLayout(view_layout)

    def setHeaderBar(self) -> None:
        header_layout = QHBoxLayout(self)
        header_layout.setAlignment(Qt.AlignmentFlag.AlignRight)

        self.searchbar = QLineEdit()
        self.searchbar.textChanged.connect(self.proxy_model.setFilterFixedString)

        self.load_tasks_button = QPushButton()
        icon_path = str(ASSETS / "refresh_icon.svg")
        self.load_tasks_button.setIcon(QIcon(icon_path))
        self.load_tasks_button.clicked.connect(self.loadTasks)

        self.settings_button = QPushButton()
        icon_path = str(ASSETS / "settings-icon-symbol.png")
        self.settings_button.setIcon(QIcon(icon_path))
        self.settings_button.clicked.connect(self.handleSettingsPressed)

        self.setSettingsActions()

        header_layout.addWidget(self.searchbar)
        header_layout.addWidget(self.load_tasks_button)
        header_layout.addWidget(self.settings_button)

        self.layout.addLayout(header_layout)

    def setSettingsActions(self) -> None:
        pass

    def updateKindsDict(self, kind: Kind) -> None:
        kind_item = QStandardItem(kind.name)
        kind_item.setData(kind, Qt.ItemDataRole.UserRole)
        self.model.appendRow(kind_item)
        self.kind_items.update({kind: kind_item})

    def addTask(self, task: Task):
        task_item = QStandardItem(task.name)
        task_item.setData(task, Qt.ItemDataRole.UserRole)

        if task.kind not in self.kind_items:
            self.updateKindsDict(task.kind)

        kind_item = self.kind_items.get(task.kind)  # get item containing corresponding kind with dictionary
        kind_item.appendRow(task_item)

    def _clearTasksList(self) -> None:
        for r in range(0, self.model.rowCount()):  # run through top level categories and remove their contents (rows)
            kind_row = self.model.item(r)
            kind_row.removeRows(0, kind_row.rowCount())

    def handleSettingsPressed(self) -> None:
        if self.settings_menu.isHidden():
            self.settings_menu.exec(self.settings_button.mapToGlobal(self.settings_button.rect().bottomLeft()))
        else:
            self.settings_menu.close()

    def handleAcceptClick(self):
        if self.view.currentIndex().row() == -1:  # no index chosen
            return

        self.handleDoubleClick()

    def handleDoubleClick(self) -> None:
        index = self.view.currentIndex()
        data = index.data(Qt.ItemDataRole.UserRole)
        if isinstance(data, Task):
            self.chooseTask(data)

    def chooseTask(self, task: Task) -> None:
        self.editor.setMode(EditorMode.SOLVE_MODE)
        self.editor.setTask(task=task)
        self.close()

    def handleDeleteClick(self) -> None:
        index = self.view.currentIndex()
        if index.row() == -1:  # no index chosen
            return

        data = index.data(Qt.ItemDataRole.UserRole)
        if not isinstance(data, Task):
            return

        self._deleteTaskFromModel(index)
        self.deleteTaskFromDatabase(data)

    def _deleteTaskFromModel(self, index: QModelIndex) -> None:
        self.model.removeRow(index.row(), index.parent())

    def deleteTaskFromDatabase(self, task: Task) -> None:
        pass

    def handleRandomTaskClick(self) -> None:
        if self.model.rowCount() == 0:  # no kinds in the model
            return

        type_id: int = randint(0, self.model.rowCount() - 1)
        type_item = self.model.itemFromIndex(self.model.index(type_id, 0))

        task_id = randint(0, type_item.rowCount() - 1)
        task_item = type_item.child(task_id, 0)
        self.chooseTask(task_item.data(Qt.ItemDataRole.UserRole))

    def loadTasks(self) -> None:
        pass
