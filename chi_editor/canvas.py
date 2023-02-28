from typing import TypeVar

from PyQt6.QtCore import QRectF
from PyQt6.QtWidgets import QGraphicsScene, QGraphicsSceneMouseEvent


Tool = TypeVar('Tool')


class Canvas(QGraphicsScene):
    current_action: Tool
    min_bounding_rect: QRectF
    """Area bounding all the items on scene"""

    def __init__(self, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs)
        self.min_bounding_rect = super().sceneRect()

    def mousePressEvent(self, event: QGraphicsSceneMouseEvent) -> None:
        self.current_action.mouse_press_event(event)

    def mouseMoveEvent(self, event: QGraphicsSceneMouseEvent) -> None:
        self.current_action.mouse_move_event(event)

    def mouseReleaseEvent(self, event: QGraphicsSceneMouseEvent) -> None:
        self.current_action.mouse_release_event(event)

    def sceneRect(self) -> QRectF:
        return self.min_bounding_rect

    def setSceneRect(self, x: float, y: float, w: float, h: float) -> None:
        self.min_bounding_rect = self.min_bounding_rect.united(QRectF(x, y, w, h))
