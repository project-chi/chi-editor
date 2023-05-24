from typing import TYPE_CHECKING

from PyQt6.QtCore import QRectF, Qt, pyqtSignal
from PyQt6.QtGui import QIcon
from PyQt6.QtWidgets import QWidget, QGraphicsView, QVBoxLayout, QPushButton, QHBoxLayout

from chi_editor.canvas import Canvas
from chi_editor.constants import ASSETS
from chi_editor.error_handler import error_handler

if TYPE_CHECKING:
    from chi_editor.toolbar.general_toolbar import GeneralToolBar


class PopupCanvas(QWidget):
    canvas: "Canvas"
    canvas_changed: pyqtSignal = pyqtSignal(str)

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

        toolbar = GeneralToolBar(canvas=self.canvas)
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
        if self.canvas.no_atoms():
            return error_handler("No molecules found fella")

        if self.canvas.more_than_one_molecule():
            return error_handler("More than one molecule you've got here fella")

        molecule = self.canvas.findMolecule()
        self.canvas_changed.emit(molecule)
