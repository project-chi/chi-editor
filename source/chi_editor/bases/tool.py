from PyQt6.QtGui import QAction, QIcon
from PyQt6.QtWidgets import QGraphicsSceneMouseEvent

from ..canvas import Canvas
from ..constants import RESOURCES


class Tool(QAction):
    def __init__(self, *args, canvas: Canvas, **kwargs) -> None:
        super().__init__(*args, **kwargs)
        self.canvas = canvas

        path = RESOURCES / 'assets' / 'toolbar' / f'{self.asset}.png'
        self.setIcon(QIcon(str(path)))
        self.setCheckable(True)
        self.triggered.connect(self.choose)

    @property
    def asset(self) -> str:
        raise NotImplemented

    def choose(self) -> None:
        self.canvas.current_action = self

    def action(self, event: QGraphicsSceneMouseEvent) -> None:
        raise NotImplemented
