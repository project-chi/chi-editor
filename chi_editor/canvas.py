from PyQt6.QtWidgets import QGraphicsScene, QGraphicsItem, QGraphicsTextItem
from PyQt6.QtGui import QFont


class Canvas(QGraphicsScene):
    shapes: list[QGraphicsItem]
    opt: str

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.setSceneRect(-100, -100, 200, 200)
        self.shapes = []
        self.opt = ""

    def setOption(self, opt):
        self.opt = opt

    def mousePressEvent(self, event):
        if self.opt == "Generate":
            text = QGraphicsTextItem()
            text.setPos(event.scenePos())
            text.setPlainText("sadasdasd")
            text.setFont(QFont("Impact"))
            text.setFlag(QGraphicsItem.GraphicsItemFlag.ItemIsMovable)
            self.shapes.append(text)
            self.addItem(text)
        else:
            super(Canvas, self).mousePressEvent(event)
