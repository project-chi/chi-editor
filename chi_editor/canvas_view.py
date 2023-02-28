from PyQt6.QtWidgets import QGraphicsView
from PyQt6.QtGui import QResizeEvent


class CanvasView(QGraphicsView):
    def resizeEvent(self, event: QResizeEvent) -> None:
        super().resizeEvent(event)
        self.scene().setSceneRect(0, 0, event.size().width(), event.size().height())
