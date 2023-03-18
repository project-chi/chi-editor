from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import QGraphicsSceneMouseEvent, QWidget, QLineEdit, QLabel, QVBoxLayout

from .structure import put_molecule
from ...bases.tool import Tool


class Smiles(Tool):
    def mouse_press_event(self, event: QGraphicsSceneMouseEvent) -> None:
        if event.button() == Qt.MouseButton.LeftButton:
            input_widget = SmilesInput()
            text, label = input_widget.initUI()

            put_molecule(self.canvas, mol, event.scenePos())

    @property
    def asset(self) -> str:
        return "structure"


class SmilesInput(QWidget):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.text_field = None
        self.saved_label = None
        self.initUI()

    def initUI(self) -> (QLineEdit, QLabel):
        self.setGeometry(300, 300, 250, 150)
        self.setWindowTitle('Interactive Text Field')

        # create a QLineEdit object and set it to the center of the window
        self.text_field = QLineEdit(self)
        self.text_field.returnPressed.connect(self.save_input)

        # create a QLabel object to display the saved text
        self.saved_label = QLabel(self)

        # create a QVBoxLayout to arrange the widgets vertically
        layout = QVBoxLayout()
        layout.addWidget(self.text_field)
        layout.addWidget(self.saved_label)
        self.setLayout(layout)

        self.show()
        return self.text_field, self.saved_label

    def save_input(self):
        text = self.text_field.text()
        self.saved_label.setText(f'Saved text: {text}')
