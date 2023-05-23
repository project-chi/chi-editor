from PyQt6.QtCore import QRectF, Qt
from PyQt6.QtGui import QIcon
from PyQt6.QtWidgets import QWidget, QGraphicsView, QVBoxLayout, QPushButton, QHBoxLayout, QMessageBox

from chi_editor.canvas import Canvas
from chi_editor.constants import ASSETS
from chi_editor.main import error_handler
from chi_editor.toolbar import CanvasToolBar


class PopupCanvas(QWidget):
    canvas: "Canvas"

    def __init__(self, *args, **kwargs) -> "None":
        super().__init__(*args, **kwargs)
        self.setWindowTitle("Write your answer here")
        self.setWindowIcon(QIcon(str(ASSETS / "project-chi.png")))
        tool_layout = QVBoxLayout(self)
        tool_layout.setContentsMargins(0, 0, 0, 0)

        view = QGraphicsView()
        view.setAlignment(Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignTop)
        view.setViewportUpdateMode(QGraphicsView.ViewportUpdateMode.FullViewportUpdate)

        self.canvas = Canvas(QRectF(100, 100, 100, 100))

        view.setScene(self.canvas)

        toolbar = CanvasToolBar(canvas=self.canvas)
        tool_layout.addWidget(toolbar)
        tool_layout.addWidget(view)

        # Button to return user input
        final_layout = QVBoxLayout()
        tool_layout.addLayout(final_layout)
        full_layout = QHBoxLayout(self)
        full_layout.setAlignment(Qt.AlignmentFlag.AlignBottom)
        return_answer = QPushButton("Return Answer", self)
        return_answer.clicked.connect(self.return_answer)
        full_layout.addWidget(return_answer)
        tool_layout.addLayout(full_layout)

    def return_answer(self) -> None:
        if self.canvas.more_than_one_molecule():
            error_handler("More than one molecule you've got here fella")
            return