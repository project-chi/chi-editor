from pathlib import Path

from PyQt6.QtWidgets import QGraphicsItem
from PyQt6.QtGui import QPixmap, QPainter, QPen, QColor


class Line(QGraphicsItem):
    texture: QPixmap
    def __init__(self, *args, vertex1: QGraphicsItem, vertex2: QGraphicsItem, **kwargs) -> None:
        pass

    def fill_texture(self):
        painter = QPainter(self.texture)
        rect = self.texture.rect()
        painter.drawLine(rect.x() + rect.width()/2, rect.y(),
                         rect.x() + rect.width()/2, rect.y() + rect.height())