from PyQt6.QtCore import Qt
from PyQt6.QtGui import QIcon
from PyQt6.QtWidgets import QMainWindow, QGraphicsView, QSizePolicy

from .canvas import Canvas
from .constants import ASSETS
from .toolbar import CanvasToolBar
from .canvas_view import CanvasView


class Editor(QMainWindow):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        canvas = Canvas(self)
        graphics_view = CanvasView(self)

        self.setWindowTitle("Project Chi")
        self.setWindowIcon(QIcon(str(ASSETS / 'project-chi.png')))
        self.resize(1000, 600)

        self.addToolBar(Qt.ToolBarArea.LeftToolBarArea, CanvasToolBar(canvas=canvas))
        self.setCentralWidget(graphics_view)

        graphics_view.setAlignment(Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignTop)
        graphics_view.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        graphics_view.setScene(canvas)
        graphics_view.setViewportUpdateMode(QGraphicsView.ViewportUpdateMode.FullViewportUpdate)
        self.show()
