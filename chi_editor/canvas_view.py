from typing import TYPE_CHECKING

from PyQt6.QtWidgets import QGraphicsView

if TYPE_CHECKING:
    from PyQt6.QtCore import QPoint
    from PyQt6.QtGui import QMouseEvent, QResizeEvent


class CanvasView(QGraphicsView):
    old_drag_pos: "QPoint"

    def resizeEvent(self, event: "QResizeEvent") -> "None":
        super().resizeEvent(event)
        self.scene().setSceneRect(
            0.0,
            0.0,
            float(event.size().width()),
            float(event.size().height())
        )

    def mousePressEvent(self, event: "QMouseEvent") -> "None":
        super().mousePressEvent(event)
        self.old_drag_pos = event.scenePosition()

    def mouseMoveEvent(self, event: "QMouseEvent") -> "None":
        if self.dragMode() != QGraphicsView.DragMode.ScrollHandDrag:
            super().mouseMoveEvent(event)
        else:
            super().mouseMoveEvent(event)
