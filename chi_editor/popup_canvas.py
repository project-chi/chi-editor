from PyQt6.QtCore import QRectF, Qt
from PyQt6.QtWidgets import QWidget, QBoxLayout, QGraphicsView, QVBoxLayout, QHBoxLayout

from chi_editor.canvas import Canvas
from chi_editor.toolbar import CanvasToolBar


class PopupCanvas(QWidget):
    def __init__(self, *args, **kwargs) -> "None":
        super().__init__(*args, **kwargs)
        tool_layout = QVBoxLayout(self)
        tool_layout.setContentsMargins(0, 0, 0, 0)

        view = QGraphicsView()
        view.setAlignment(Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignTop)
        view.setViewportUpdateMode(QGraphicsView.ViewportUpdateMode.FullViewportUpdate)
        tool_layout.addWidget(view)

        canvas = Canvas(QRectF(100, 100, 100, 100))

        view.setScene(canvas)

        toolbar = CanvasToolBar(canvas=canvas)
        tool_layout.addWidget(toolbar)

        final_layout = QVBoxLayout()
        tool_layout.addLayout(final_layout)


