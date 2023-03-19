from enum import Enum

from PyQt6.QtCore import Qt, QRectF
from PyQt6.QtGui import QIcon, QTransform
from PyQt6.QtWidgets import QMainWindow, QGraphicsView, QPushButton, QVBoxLayout, QStackedWidget, QToolBar, QWidget

from .canvas import Canvas
from .constants import ASSETS
from .toolbar import CanvasToolBar
from .menubar.menubar import CanvasMenuBar

from .editor_mode import EditorMode


class Editor(QMainWindow):
    # Hierarchy:
    #
    #   window:
    #       workspace:
    #           graphics_view:
    #               canvas:
    #       toolbar:

    # Contains everything except toolbar
    workspace: QStackedWidget

    # QGraphicsView contains drawable space
    graphics_view: QGraphicsView

    # GraphicsScene where to draw all graphical objects
    canvas: Canvas

    # Toolbar that contains tools for manipulating canvas in free mode
    _drawing_tool_bar: QToolBar

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        # Window settings
        self.setWindowTitle("Project Chi")
        self.setWindowIcon(QIcon(str(ASSETS / 'project-chi.png')))
        self.resize(1000, 600)

        # The biggest part of interface
        self.workspace = QStackedWidget()  # create workspace
        self.setCentralWidget(self.workspace)

        # Set default (free) mode
        self.createModes()
        self.workspace.widget(0).show()

        # Add custom menuBar
        menubar = CanvasMenuBar(editor=self)
        self.setMenuBar(menubar)

        # Add toolbar
        self._drawing_tool_bar = CanvasToolBar(canvas=self.canvas, parent=self)
        self.addToolBar(Qt.ToolBarArea.LeftToolBarArea, self._drawing_tool_bar)

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

    def setMode(self, mode: EditorMode) -> None:
        self.workspace.widget(self.workspace.currentIndex()).hide()
        self.workspace.widget(mode.value).show()

    def createModes(self) -> None:
        self.workspace.addWidget(self.getFreeModeLayout())
        self.workspace.addWidget(self.getSolveModeLayout())
        self.workspace.addWidget(self.getCreateModeLayout())

    def getFreeModeLayout(self) -> QWidget:
        # Initialize QGraphicsView
        self.graphics_view = QGraphicsView(self)  # create QGraphicsView
        self.graphics_view \
            .setAlignment(Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignTop)  # Set QGraphicsView position
        self.graphics_view.setViewportUpdateMode(QGraphicsView.ViewportUpdateMode.FullViewportUpdate)

        # Initialize GGraphicsScene called canvas
        self.canvas = Canvas(QRectF(self.graphics_view.geometry()))

        # Bind GraphicsScene to GraphicsView
        self.graphics_view.setScene(self.canvas)

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

        self.graphics_view.setLayout(v_button_group)

        free_mode_widget = QWidget()
        layout = QVBoxLayout(free_mode_widget)
        layout.addWidget(self.graphics_view)
        return free_mode_widget

    def getSolveModeLayout(self) -> QWidget:
        solve_mode_widget = QWidget()
        layout = QVBoxLayout(solve_mode_widget)
        return solve_mode_widget

    def getCreateModeLayout(self) -> QWidget:
        create_mode_widget = QWidget()
        layout = QVBoxLayout(create_mode_widget)
        return create_mode_widget
