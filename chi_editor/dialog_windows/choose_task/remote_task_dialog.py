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
    # Tasks database server
    server: Server

    def __init__(self, *args, editor: "Editor", **kwargs) -> None:
        super().__init__(*args, editor=editor, **kwargs)
        self.server = Server(default_url)

    def loadTasks(self) -> None:
        self._clearTasksList()

        tasks = self.server.get_tasks_raw()

        for task in tasks:
            self.addTask(task)
