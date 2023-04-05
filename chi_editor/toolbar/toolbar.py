from typing import TYPE_CHECKING

from PyQt6.QtCore import QSize
from PyQt6.QtGui import QActionGroup
from PyQt6.QtWidgets import QToolBar

from chi_editor.toolbar.tools import tools

if TYPE_CHECKING:
    from PyQt6.QtGui import QAction

    from chi_editor.canvas import Canvas


class CanvasToolBar(QToolBar):
    _canvas: "Canvas"

    def __init__(self, *args, canvas: "Canvas", **kwargs) -> "None":
        super().__init__(*args, **kwargs)
        self._canvas = canvas
        self.setStyleSheet("""QToolBar { background-color: rgb(212, 204, 234); }""")
        self.setMovable(False)

        action_group = QActionGroup(self)
        action_group.setExclusive(True)

        for Tool in tools:
            tool = Tool(canvas)
            self.addAction(tool)
            action_group.addAction(tool)
            self.widgetForAction(tool).setStyleSheet("padding: 5px")

        self.actionTriggered.connect(self.changeAction)

    def changeAction(self, action: "QAction") -> "None":
        self._canvas.current_action = action
