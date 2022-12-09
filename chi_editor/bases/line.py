from pathlib import Path

from PyQt6.QtWidgets import QGraphicsItem


class Line(QGraphicsItem):
    def __init__(self, *args, vertex1: QGraphicsItem, vertex2: QGraphicsItem, **kwargs) -> None:
        pass