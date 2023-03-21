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
    view_free: QGraphicsView
    view_solver: QGraphicsView

    # GraphicsScene where to draw all graphical objects
    canvas_free: Canvas
    canvas_solver: Canvas

    # Toolbar that contains tools for manipulating canvas in free mode
    toolbars: list[QToolBar]
    mode: EditorMode

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
        _drawing_tool_bar = CanvasToolBar(canvas=self.canvas_free, parent=self)
        _solver_tool_bar = CanvasToolBar(canvas=self.canvas_solver, parent=self)
        self.toolbars = [_drawing_tool_bar, _solver_tool_bar, None]
        self.addToolBar(Qt.ToolBarArea.LeftToolBarArea, _drawing_tool_bar)
        self.addToolBar(Qt.ToolBarArea.LeftToolBarArea, _solver_tool_bar)
        for toolbar in self.toolbars:
            if toolbar:
                toolbar.toggleViewAction().setChecked(False)
                toolbar.toggleViewAction().trigger()
        _solver_tool_bar.toggleViewAction().trigger()
        self.mode = 0

    def zoom_in(self):
        # Get the current scale factor of the view
        current_scale = self.view_free.transform().m11()

        # Update the scale factor of the view
        new_scale = current_scale * 1.2
        self.view_free.setTransform(QTransform.fromScale(new_scale, new_scale))

    def zoom_out(self):
        # Get the current scale factor of the view
        current_scale = self.view_free.transform().m11()

        # Update the scale factor of the view
        new_scale = current_scale / 1.2
        self.view_free.setTransform(QTransform.fromScale(new_scale, new_scale))

    def setMode(self, mode: EditorMode) -> None:
        self.toolbars[self.mode].toggleViewAction().trigger()
        self.workspace.widget(self.workspace.currentIndex()).hide()
        self.workspace.widget(mode.value).show()
        self.toolbars[mode.value].toggleViewAction().trigger()
        self.mode = mode.value

    def createModes(self) -> None:
        self.workspace.addWidget(self.getFreeModeLayout())
        self.workspace.addWidget(self.getSolveModeLayout())
        self.workspace.addWidget(self.getCreateModeLayout())

    def getFreeModeLayout(self) -> QWidget:
        # Initialize QGraphicsView
        self.view_free = QGraphicsView(self)  # create QGraphicsView
        self.view_free \
            .setAlignment(Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignTop)  # Set QGraphicsView position
        self.view_free.setViewportUpdateMode(QGraphicsView.ViewportUpdateMode.FullViewportUpdate)

        # Initialize GGraphicsScene called canvas
        self.canvas_free = Canvas(QRectF(self.view_free.geometry()))

        # Bind GraphicsScene to GraphicsView
        self.view_free.setScene(self.canvas_free)

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

        self.view_free.setLayout(v_button_group)

        free_mode_widget = QWidget()
        layout = QVBoxLayout(free_mode_widget)
        layout.addWidget(self.view_free)
        return free_mode_widget

    def getSolveModeLayout(self) -> QWidget:
        # Initialize QGraphicsView
        self.view_solver = QGraphicsView(self)  # create QGraphicsView
        self.view_solver \
            .setAlignment(Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignTop)  # Set QGraphicsView position
        self.view_solver.setViewportUpdateMode(QGraphicsView.ViewportUpdateMode.FullViewportUpdate)

        # Initialize GGraphicsScene called canvas
        self.canvas_solver = Canvas(QRectF(self.view_solver.geometry()))

        # Bind GraphicsScene to GraphicsView
        self.view_solver.setScene(self.canvas_solver)

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

        self.view_solver.setLayout(v_button_group)

        solver_mode_widget = QWidget()
        layout = QVBoxLayout(solver_mode_widget)
        layout.addWidget(self.view_solver)
        return solver_mode_widget

    def getCreateModeLayout(self) -> QWidget:
        create_mode_widget = QWidget()
        layout = QVBoxLayout(create_mode_widget)
        return create_mode_widget
