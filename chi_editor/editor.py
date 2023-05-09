from PyQt6.QtCore import QRectF, Qt
from PyQt6.QtGui import QIcon, QTransform
from PyQt6.QtWidgets import (
    QGraphicsView,
    QHBoxLayout,
    QMainWindow,
    QPushButton,
    QVBoxLayout,
    QWidget,
    QStackedWidget,
    QToolBar,
    QSizePolicy,
    QLabel,
    QGroupBox,
)

from rdkit import Chem

from chi_editor.canvas import Canvas
from chi_editor.constants import ASSETS
from chi_editor.task_creator import InputDialog
from chi_editor.toolbar import CanvasToolBar
from chi_editor.menubar.menubar import CanvasMenuBar
from chi_editor.dialog_windows.task_result_dialog import TaskResultDialog

from chi_editor.api.task import Task

from chi_editor.editor_mode import EditorMode


class Editor(QMainWindow):
    # Hierarchy:
    #
    #   window:
    #       workspace:
    #           graphics_view:
    #               canvas:
    #       toolbar:

    # Contains everything except toolbar
    workspace: "QStackedWidget"

    # QGraphicsView contains drawable space
    views: list[QGraphicsView] = [None, None, None]

    # GraphicsScene where to draw all graphical objects
    canvases: list[Canvas] = [None, None, None]

    # Toolbar that contains tools for manipulating canvas in free mode
    toolbars: list[QToolBar] = [None, None, None]

    # Dialog with solving results
    result_dialog: TaskResultDialog

    # Task that user currently solves
    task: Task = None

    # Text of the formulation of the current task
    formulation: QLabel = None

    def __init__(self, *args, **kwargs) -> "None":
        super().__init__(*args, **kwargs)

        # Window settings
        self.setWindowTitle("Project Chi")
        self.setWindowIcon(QIcon(str(ASSETS / "project-chi.png")))
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

        # Add toolbars
        for index in range(3):
            self.toolbars[index] = CanvasToolBar(canvas=self.canvases[index], parent=self)
            self.addToolBar(Qt.ToolBarArea.LeftToolBarArea, self.toolbars[index])
            self.zoom_out(index)
            self.zoom_out(index)
        for toolbar in self.toolbars:
            toolbar.toggleViewAction().setChecked(False)
            toolbar.toggleViewAction().trigger()
            toolbar.toggleViewAction().trigger()
        self.toolbars[0].toggleViewAction().trigger()

        # Create dialogs (they won't show instantly)

        self.result_dialog = TaskResultDialog(editor=self)

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

    def setTask(self, task: Task) -> None:
        self.task = task
        self.setFormulationOfTask()

    def createModes(self) -> None:
        self.workspace.addWidget(self.getFreeModeLayout())
        self.workspace.addWidget(self.getSolveModeLayout())
        self.workspace.addWidget(self.getCreateModeLayout())

    def getFreeModeLayout(self) -> QWidget:
        return self.getLayout(0)

    def getSolveModeLayout(self) -> QWidget:
        solve_widget = self.getLayout(1)
        solve_layout = solve_widget.layout()

        # Button to submit answers
        submit = QPushButton("Submit", self)
        submit.setFixedSize(submit.sizeHint())
        submit.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed)
        submit.clicked.connect(self.submitAnswer)

        # Text box creation

        box_layout = QVBoxLayout()
        group_box = QGroupBox("Formulation")
        formulation = QLabel("The formulation of a task will appear here")
        self.formulation = formulation
        box_layout.addWidget(formulation)
        box_layout.setAlignment(Qt.AlignmentFlag.AlignTop)
        group_box.setLayout(box_layout)

        # Formulation of the task
        solve_layout.setAlignment(Qt.AlignmentFlag.AlignCenter)
        solve_layout.addWidget(submit, alignment=Qt.AlignmentFlag.AlignBottom)
        solve_layout.addWidget(group_box, alignment=Qt.AlignmentFlag.AlignTop)

        return solve_widget

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

        # Main solve canvas layout
        solve_canvas_layout = QHBoxLayout(self)

        # Box contains magnifying glass
        zoom_layout = QVBoxLayout(self)
        zoom_layout.setAlignment(Qt.AlignmentFlag.AlignRight)

        scale_plus = QPushButton("Zoom In", self)
        scale_plus.clicked.connect(lambda: self.zoom_in(index))

        scale_minus = QPushButton("Zoom Out", self)
        scale_minus.clicked.connect(lambda: self.zoom_out(index))

        zoom_layout.addWidget(scale_minus)
        zoom_layout.addWidget(scale_plus)

        # Stack layouts
        solve_canvas_layout.addLayout(zoom_layout)

        self.views[index].setLayout(solve_canvas_layout)

        solver_mode_widget = QWidget()
        layout = QVBoxLayout(solver_mode_widget)
        layout.addWidget(self.views[index])
        return solver_mode_widget

    def changeCanvas(self, canvas: Canvas, index: int):
        self.views[index].setScene(canvas)

    def getCreateModeLayout(self) -> QWidget:
        create_widget = self.getLayout(EditorMode.CREATE_MODE.value)
        create_layout = create_widget.layout()
        create_mode_widget = InputDialog(self.canvases[EditorMode.CREATE_MODE.value])
        QVBoxLayout(create_mode_widget)
        create_layout.addWidget(create_mode_widget)
        return create_widget

    def submitAnswer(self) -> None:
        smiles_answer = self.canvases[EditorMode.SOLVE_MODE.value].findMolecule()
        answer_is_correct = self.checkAnswer(smiles_answer)

        if answer_is_correct:
            self.openResultDialog("Correct")
        else:
            self.openResultDialog("Wrong")

    def openResultDialog(self, message: str) -> None:
        self.result_dialog.setText(message)
        self.result_dialog.exec()

    def setFormulationOfTask(self) -> None:
        self.formulation.setText(self.task.problem)

    def checkAnswer(self, user_answer: str) -> bool:
        mol_user = Chem.MolFromSmiles(user_answer)
        if mol_user == None:
            return False
        Chem.RemoveStereochemistry(mol_user)
        user_smiles = Chem.MolToSmiles(mol_user)
        mol_correct = Chem.MolFromSmiles(self.task.solution)
        Chem.RemoveStereochemistry(mol_correct)
        correct_smiles = Chem.MolToSmiles(mol_correct)
        canon_user_answer = Chem.CanonSmiles(user_smiles)
        correct_answer = Chem.CanonSmiles(correct_smiles)
        return correct_answer == canon_user_answer
