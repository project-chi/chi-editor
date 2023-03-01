from PyQt6.QtCore import Qt, QRectF
from PyQt6.QtGui import QIcon, QColor, QTransform, QBrush
from PyQt6.QtWidgets import QMainWindow, QGraphicsView, QSizePolicy, QPushButton, QVBoxLayout, QHBoxLayout, QWidget, \
    QStackedLayout, QBoxLayout

from .canvas import Canvas
from .constants import ASSETS
from .toolbar import CanvasToolBar
from .canvas_view import CanvasView


class Editor(QMainWindow):

    # Hierarchy:
    #
    #   window:
    #       workspace:
    #           graphics_view:
    #               canvas:
    #       toolbar:

    # Contains everything except toolbar
    workspace: QWidget

    # QGraphicsView contains drawable space
    graphics_view: QGraphicsView

    # GraphicsScene where to draw all graphical objects
    canvas: Canvas

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        # Window settings
        self.setWindowTitle("Project Chi")
        self.setWindowIcon(QIcon(str(ASSETS / 'project-chi.png')))
        self.resize(1000, 600)

        # The biggest part of interface
        self.workspace = QWidget()  # create workspace

        # Initialize QGraphicsView
        self.graphics_view = QGraphicsView(self)    # create QGraphicsView
        self.graphics_view\
            .setAlignment(Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignTop)     # Set QGraphicsView position

        # Initialize GGraphicsScene called canvas
        self.canvas = Canvas(self)
        self.canvas.setSceneRect(QRectF(self.geometry()))

        # Bind GraphicsScene to GraphicsView
        self.graphics_view.setScene(self.canvas)

        h_button_group = QHBoxLayout(self)

        # Box contains magnifying glass
        v_button_group = QVBoxLayout(self)
        v_button_group.addStretch(1)
        v_button_group.setAlignment(Qt.AlignmentFlag.AlignRight)

        scale_plus = QPushButton("Zoom In", self)
        scale_plus.clicked.connect(self.zoom_in)

        scale_minus = QPushButton("Zoom Out", self)
        scale_minus.clicked.connect(self.zoom_out)

        v_button_group.addWidget(scale_minus)
        v_button_group.addWidget(scale_plus)
        v_button_group.addStretch(1)

        h_button_group.addStretch()
        h_button_group.addLayout(v_button_group)
        h_button_group.setAlignment(Qt.AlignmentFlag.AlignCenter)

        layout = QVBoxLayout(self)
        layout.addWidget(self.graphics_view)

        self.graphics_view.setLayout(h_button_group)

        self.workspace.setLayout(layout)
        self.setCentralWidget(self.workspace)

        # Add left toolbar
        self.addToolBar(Qt.ToolBarArea.LeftToolBarArea, CanvasToolBar(canvas=self.canvas))

    def zoom_in(self):
        # Get the current scale factor of the view
        current_scale = self.graphics_view.transform().m11()

        # Update the scale factor of the view
        new_scale = current_scale * 1.2
        self.graphics_view.setTransform(QTransform.fromScale(new_scale, new_scale))

    def zoom_out(self):
        # Get the current scale factor of the view
        current_scale = self.graphics_view.transform().m11()

        # Update the scale factor of the view
        new_scale = current_scale / 1.2
        self.graphics_view.setTransform(QTransform.fromScale(new_scale, new_scale))
