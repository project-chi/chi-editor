from enum import Enum

from PyQt6.QtCore import Qt, QRectF
from PyQt6.QtGui import QIcon, QTransform
from PyQt6.QtWidgets import QMainWindow, QGraphicsView, QPushButton, QVBoxLayout, QStackedWidget, QToolBar, QWidget, \
    QHBoxLayout, QSizePolicy

from rdkit import Chem

from .canvas import Canvas
from .constants import ASSETS
from .toolbar import CanvasToolBar
from .menubar.menubar import CanvasMenuBar
from .choose_task_dialog import ChooseTaskDialog
from .task_result_dialog import TaskResultDialog

from .tasks.task import Task
from .bases.molecule.molecule import Molecule
from .bases.alpha_atom import AlphaAtom
from .playground import mol_from_graphs

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

    # Dialog where user chooses a task to solve
    choose_task_dialog: ChooseTaskDialog

    # Dialog with solving results
    result_dialog: TaskResultDialog

    # Task that user currently solves
    task: Task = None

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

        # Create dialogs (they won't show for now)
        self.choose_task_dialog = ChooseTaskDialog(editor=self)
        self.choose_task_dialog.setModal(True)

        self.result_dialog = TaskResultDialog(editor=self)
        self.result_dialog.setModal(True)

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

    def setTask(self, task: Task) -> None:
        self.task = task

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

        # Main solve canvas layout
        solve_canvas_layout = QHBoxLayout(self)

        # Box contains magnifying glass
        zoom_layout = QVBoxLayout(self)
        zoom_layout.setAlignment(Qt.AlignmentFlag.AlignRight)

        scale_plus = QPushButton("Zoom In", self)
        scale_plus.clicked.connect(self.zoom_in)

        scale_minus = QPushButton("Zoom Out", self)
        scale_minus.clicked.connect(self.zoom_out)

        zoom_layout.addWidget(scale_minus)
        zoom_layout.addWidget(scale_plus)

        # Button to submit answers
        submit = QPushButton("Submit", self)
        submit.setFixedSize(submit.sizeHint())
        submit.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed)
        submit.clicked.connect(self.submitAnswer)

        # Stack layouts
        solve_canvas_layout.addWidget(submit)
        solve_canvas_layout.addLayout(zoom_layout)

        self.view_solver.setLayout(solve_canvas_layout)

        solver_mode_widget = QWidget()
        layout = QVBoxLayout(solver_mode_widget)
        layout.addWidget(self.view_solver)
        return solver_mode_widget

    def submitAnswer(self) -> None:
        items = self.canvas_solver.items()

        if len(items) == 0:
            self.openResultDialog("No answer found")
            return

        molecule: Molecule = None

        # Простите пожалуйста, Please forgive me, Entschuldigung bitte
        for item in items:
            if isinstance(item, AlphaAtom):
                molecule = item.molecule
                break

        if molecule is None:
            self.openResultDialog("Something was found, but not a molecule")
            return

        smiles_answer = Chem.MolToSmiles(mol_from_graphs(molecule))
        answer_is_correct = self.task.checkAnswer(smiles_answer)

        if answer_is_correct:
            self.openResultDialog("Correct")
        else:
            self.openResultDialog("Wrong")

    def getCreateModeLayout(self) -> QWidget:
        create_mode_widget = QWidget()
        layout = QVBoxLayout(create_mode_widget)
        return create_mode_widget

    def openChooseTaskDialog(self) -> None:
        self.choose_task_dialog.exec()

    def openResultDialog(self, message: str) -> None:
        self.result_dialog.setText(message)
        self.result_dialog.exec()
