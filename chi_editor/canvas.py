from PyQt6.QtWidgets import QGraphicsScene, QGraphicsItem, QGraphicsTextItem
from PyQt6.QtGui import QFont, QAction


class Canvas(QGraphicsScene):
    shapes: list[QGraphicsItem]
    opt: str
    mouse_press_event: None

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # self.setSceneRect(-100, -100, 200, 200)
        self.shapes = []
        self.opt = ""
        self.mouse_press_event = self.textEvent

    def setText(self):
        self.mouse_press_event = self.textEvent

    def setArrow(self):
        self.mouse_press_event = self.arrowEvent

    def mousePressEvent(self, event):
        self.mouse_press_event(event)
        # if self.opt == "Text":
        #     text = QGraphicsTextItem()
        #     text.setPos(event.scenePos())
        #     text.setPlainText("sadasdasd")
        #     text.setFont(QFont("Impact"))
        #     text.setFlag(QGraphicsItem.GraphicsItemFlag.ItemIsMovable)
        #     self.shapes.append(text)
        #     self.addItem(text)
        # else:
        #     super(Canvas, self).mousePressEvent(event)

    def arrowEvent(self, event):
        super(Canvas, self).mousePressEvent(event)

    def textEvent(self, event):
        text = QGraphicsTextItem()
        text.setPos(event.scenePos())
        text.setPlainText("sadasdasd")
        text.setFont(QFont("Impact"))
        text.setFlag(QGraphicsItem.GraphicsItemFlag.ItemIsMovable)
        self.shapes.append(text)
        self.addItem(text)
