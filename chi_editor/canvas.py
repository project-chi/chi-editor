from typing import TypeVar

from PyQt6.QtCore import QRectF
from PyQt6.QtWidgets import QGraphicsScene, QGraphicsSceneMouseEvent


Tool = TypeVar('Tool')


class Canvas(QGraphicsScene):
    current_action: Tool

    def __init__(self, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs)

    def mousePressEvent(self, event: QGraphicsSceneMouseEvent) -> None:
        self.current_action.mouse_press_event(event)

    def mouseMoveEvent(self, event: QGraphicsSceneMouseEvent) -> None:
        self.current_action.mouse_move_event(event)

    def mouseReleaseEvent(self, event: QGraphicsSceneMouseEvent) -> None:
        self.current_action.mouse_release_event(event)

    def minSceneRect(self) -> QRectF:
        return super().sceneRect()
