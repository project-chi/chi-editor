from typing import TypeVar
from typing import overload

from PyQt6.QtCore import QRectF
from PyQt6.QtWidgets import QGraphicsScene, QGraphicsSceneMouseEvent


Tool = TypeVar('Tool')


class Canvas(QGraphicsScene):
    current_action: Tool
    min_scene_rect: QRectF

    def __init__(self, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs)
        self.min_scene_rect = super().sceneRect()

    def mousePressEvent(self, event: QGraphicsSceneMouseEvent) -> None:
        self.current_action.mouse_press_event(event)

    def mouseMoveEvent(self, event: QGraphicsSceneMouseEvent) -> None:
        self.current_action.mouse_move_event(event)

    def mouseReleaseEvent(self, event: QGraphicsSceneMouseEvent) -> None:
        self.current_action.mouse_release_event(event)

    def sceneRect(self) -> QRectF:
        return self.min_scene_rect

    @overload
    def setSceneRect(self, sceneRect: QRectF) -> None:
        ...

    def enlargeScene(self, sceneRect: QRectF) -> None:
        self.min_scene_rect = self.min_scene_rect.united(sceneRect)

    @overload
    def setSceneRect(self, x: float, y: float, w: float, h: float) -> None:
        ...

    def setSceneRect(self, *args) -> None:
        match args:
            case [QRectF() as sceneRect]:
                self.enlargeScene(sceneRect)
            case [float() as x, float() as y, float() as w, float() as h]:
                self.enlargeScene(QRectF(x, y, w, h))
            case _:
                raise TypeError("wrong signature")