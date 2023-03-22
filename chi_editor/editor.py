from enum import Enum

from PyQt6.QtCore import Qt, QRectF
from PyQt6.QtGui import QIcon, QTransform
from PyQt6.QtWidgets import QMainWindow, QGraphicsView, QPushButton, QVBoxLayout, QStackedWidget, QToolBar, QWidget

from .canvas import Canvas
from .canvas_view import CanvasView
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
    views: list[QGraphicsView] = [None, None, None]

    # GraphicsScene where to draw all graphical objects
    canvases: list[Canvas] = [None, None, None]

    # Toolbar that contains tools for manipulating canvas in free mode
    toolbars: list[QToolBar] = [None, None, None]

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
        for index in range(3):
            self.toolbars[index] = CanvasToolBar(canvas=self.canvases[index], parent=self)
            self.addToolBar(Qt.ToolBarArea.LeftToolBarArea, self.toolbars[index])
        for toolbar in self.toolbars:
            toolbar.toggleViewAction().setChecked(False)
            toolbar.toggleViewAction().trigger()
            toolbar.toggleViewAction().trigger()
        self.toolbars[0].toggleViewAction().trigger()

    def zoom_in(self, index: int):
        # Get the current scale factor of the view
        current_scale = self.views[index].transform().m11()

        # Update the scale factor of the view
        new_scale = current_scale * 1.2
        self.views[index].setTransform(QTransform.fromScale(new_scale, new_scale))

    def zoom_out(self, index: int):
        # Get the current scale factor of the view
        current_scale = self.views[index].transform().m11()

        # Update the scale factor of the view
        new_scale = current_scale / 1.2
        self.views[index].setTransform(QTransform.fromScale(new_scale, new_scale))

    def setMode(self, mode: EditorMode) -> None:
        self.toolbars[self.workspace.currentIndex()].toggleViewAction().trigger()
        self.workspace.setCurrentIndex(mode.value)
        self.toolbars[self.workspace.currentIndex()].toggleViewAction().trigger()

    def createModes(self) -> None:
        self.workspace.addWidget(self.getFreeModeLayout())
        self.workspace.addWidget(self.getSolveModeLayout())
        self.workspace.addWidget(self.getCreateModeLayout())

    def getFreeModeLayout(self) -> QWidget:
        return self.getLayout(0)

    def getSolveModeLayout(self) -> QWidget:
        return self.getLayout(1)

    def getLayout(self, index) -> QWidget:
        # Initialize QGraphicsView
        self.views[index] = QGraphicsView(self)  # create QGraphicsView
        self.views[index] \
            .setAlignment(Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignTop)  # Set QGraphicsView position
        self.views[index].setViewportUpdateMode(QGraphicsView.ViewportUpdateMode.FullViewportUpdate)

        # Initialize GGraphicsScene called canvas
        self.canvases[index] = Canvas(QRectF(self.views[index].geometry()))

        # Bind GraphicsScene to GraphicsView
        self.views[index].setScene(self.canvases[index])

        # Box contains magnifying glass
        v_button_group = QVBoxLayout(self)
        v_button_group.addStretch(1)
        v_button_group.setAlignment(Qt.AlignmentFlag.AlignRight)

        scale_plus = QPushButton("Zoom In", self)
        scale_plus.clicked.connect(lambda: self.zoom_in(index))

        scale_minus = QPushButton("Zoom Out", self)
        scale_minus.clicked.connect(lambda: self.zoom_out(index))

        v_button_group.addWidget(scale_minus)
        v_button_group.addWidget(scale_plus)
        v_button_group.addStretch(1)

        self.views[index].setLayout(v_button_group)

        solver_mode_widget = QWidget()
        layout = QVBoxLayout(solver_mode_widget)
        layout.addWidget(self.views[index])
        return solver_mode_widget

    def getCreateModeLayout(self) -> QWidget:
        create_mode_widget = QWidget()
        layout = QVBoxLayout(create_mode_widget)
        return create_mode_widget
