from PyQt6.QtCore import Qt, QRectF
from PyQt6.QtGui import QIcon
from PyQt6.QtWidgets import QMainWindow, QGraphicsView, QSizePolicy

from .canvas import Canvas
from .constants import ASSETS
from .toolbar import CanvasToolBar


class Editor(QMainWindow):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        canvas = Canvas(self)
        canvas.setSceneRect(QRectF(self.geometry()))

        graphics_view = QGraphicsView(self)
        graphics_view.setAlignment(Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignTop)
        graphics_view.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        graphics_view.setScene(canvas)
        graphics_view.setMaximumSize(500, 300)
        graphics_view.setViewportUpdateMode(QGraphicsView.ViewportUpdateMode.FullViewportUpdate)

        self.setWindowTitle("Project Chi")
        self.setWindowIcon(QIcon(str(ASSETS / 'project-chi.png')))
        self.resize(1000, 800)

        self.addToolBar(Qt.ToolBarArea.LeftToolBarArea, CanvasToolBar(canvas=canvas))
        self.setCentralWidget(graphics_view)
