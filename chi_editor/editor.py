from PyQt6.QtCore import Qt, QRectF
from PyQt6.QtGui import QIcon, QColor, QTransform
from PyQt6.QtWidgets import QMainWindow, QGraphicsView, QSizePolicy, QPushButton, QVBoxLayout, QHBoxLayout, QWidget, \
    QStackedLayout

from .canvas import Canvas
from .constants import ASSETS
from .toolbar import CanvasToolBar


class Editor(QMainWindow):
    canvas: Canvas
    graphics_view: QGraphicsView

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.canvas = Canvas(self)
        self.canvas.setSceneRect(QRectF(self.geometry()))

        self.graphics_view = QGraphicsView(self)
        self.graphics_view.setAlignment(Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignTop)
        self.graphics_view.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        self.graphics_view.setScene(self.canvas)
        self.graphics_view.setViewportUpdateMode(QGraphicsView.ViewportUpdateMode.FullViewportUpdate)

        self.setWindowTitle("Project Chi")
        self.setWindowIcon(QIcon(str(ASSETS / 'project-chi.png')))
        self.resize(1000, 800)

        self.addToolBar(Qt.ToolBarArea.LeftToolBarArea, CanvasToolBar(canvas=self.canvas))
        # self.setCentralWidget(self.graphics_view)

        # Create a button that changes the scale of the view
        scale_plus = QPushButton("Zoom In", self)
        scale_plus.clicked.connect(self.zoom_in)

        scale_minus = QPushButton("Zoom Out", self)
        scale_minus.clicked.connect(self.zoom_out)

        # Create a vertical box layout to hold the buttons
        vbox = QVBoxLayout()
        vbox.addStretch()
        vbox.addWidget(scale_minus)
        vbox.addWidget(scale_plus)

        # Create a horizontal box layout to hold the vertical layout
        hbox = QHBoxLayout()
        hbox.addStretch()
        hbox.addLayout(vbox)

        layout = QVBoxLayout()
        layout.addWidget(self.graphics_view)
        layout.addLayout(hbox)

        # Create a central widget and set its layout
        central_widget = QWidget()
        central_widget.setLayout(layout)

        # Set the central widget of the QMainWindow
        self.setCentralWidget(central_widget)

        # Set the view to be interactable
        self.graphics_view.setInteractive(True)

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