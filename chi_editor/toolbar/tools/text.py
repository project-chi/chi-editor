from PyQt6.QtGui import QFont, QFocusEvent
from PyQt6.QtWidgets import QGraphicsSceneMouseEvent, QGraphicsTextItem, QGraphicsItem, QLineEdit

from ...bases.tool import Tool


class MyLineEdit(QLineEdit):
    def focusOutEvent(self, a0: QFocusEvent) -> None:
        self.close()


class Text(Tool):
    def mouse_press_event(self, event: QGraphicsSceneMouseEvent) -> None:
        text_line = MyLineEdit()
        text_line.show()
        text_line.setFocus()
        text_line.move(event.pos().toPoint())
        text_line.returnPressed.connect(lambda: self.read_line(text_line))
        # new_text = QGraphicsTextItem()
        # new_text.setPos(event.scenePos())
        # new_text.setPlainText('H')
        # new_text.setFont(QFont('Impact'))
        # new_text.setFlag(QGraphicsItem.GraphicsItemFlag.ItemIsMovable)
        #
        # self.canvas.addItem(new_text)

    def read_line(self, text_line: QLineEdit) -> None:
        new_text = QGraphicsTextItem()
        new_text.setPos(text_line.pos().toPointF())
        new_text.setPlainText(text_line.text())
        new_text.setFont(QFont('Impact'))
        new_text.setFlag(QGraphicsItem.GraphicsItemFlag.ItemIsMovable)

        self.canvas.addItem(new_text)
        text_line.close()

    @property
    def asset(self) -> str:
        return 'text'

