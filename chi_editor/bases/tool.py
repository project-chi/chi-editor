from pathlib import Path

from PyQt6.QtGui import QAction, QIcon
from PyQt6.QtWidgets import QGraphicsSceneMouseEvent

from ..canvas import Canvas
from ..constants import RESOURCES


class Tool(QAction):
    def __init__(self, *args, canvas: Canvas, **kwargs) -> None:
        super().__init__(*args, **kwargs)
        self.canvas = canvas

        icon_path = RESOURCES / 'assets' / 'toolbar' / f'{self.asset}.png'
        if not icon_path.exists():
            raise ValueError(
                f'can\'t find assets for {self.__class__.__name__}, checked path is {icon_path}\n'
                f'    module reference: {self.__module__}'
            )

        self.setIcon(QIcon(str(icon_path)))
        self.setCheckable(True)
        self.triggered.connect(self.choose)

    def choose(self) -> None:
        self.canvas.current_action = self

    def mouse_press_event(self, event: QGraphicsSceneMouseEvent) -> None:
        pass

    def mouse_move_event(self, event: QGraphicsSceneMouseEvent) -> None:
        pass

    def mouse_release_event(self, event: QGraphicsSceneMouseEvent) -> None:
        pass

    @property
    def asset(self) -> str:
        """Returns the name of the icon for this tool without extension."""
        raise NotImplemented
