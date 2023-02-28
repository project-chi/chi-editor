from PyQt6.QtWidgets import QGraphicsView
from PyQt6.QtGui import QResizeEvent


class CanvasView(QGraphicsView):
    def resizeEvent(self, event: QResizeEvent) -> None:
        super().resizeEvent(event)
        self.scene().setSceneRect(float(0), float(0), float(event.size().width()), float(event.size().height()))
