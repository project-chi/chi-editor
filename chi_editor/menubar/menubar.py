from PyQt6.QtGui import QAction, QActionGroup
from PyQt6.QtWidgets import QMenuBar

from .menu_tools import menu_tools
from ..canvas import Canvas


class CanvasMenuBar(QMenuBar):
    _canvas: Canvas

    def __init__(self, *args, canvas: Canvas, **kwargs) -> None:
        super().__init__(*args, **kwargs)
        self._canvas = canvas
