from PyQt6.QtGui import QAction, QActionGroup
from PyQt6.QtWidgets import QToolBar, QGraphicsView

from ..toolbar.tools.drag import Drag
from .tools import tools
from ..canvas import Canvas


class CanvasToolBar(QToolBar):
    _canvas: Canvas

    def __init__(self, *args, canvas: Canvas, **kwargs) -> None:
        super().__init__(*args, **kwargs)
        self._canvas = canvas
        self.setStyleSheet("""QToolBar { background-color: rgb(212, 204, 234); }""")
        self.setMovable(False)

        action_group = QActionGroup(self)
        action_group.setExclusive(True)

        for Tool in tools:
            tool = Tool(canvas=canvas)
            self.addAction(tool)
            action_group.addAction(tool)

        self.actionTriggered.connect(self.changeAction)

    def changeAction(self, action: QAction) -> None:
        self._canvas.current_action = action