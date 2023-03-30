from PyQt6.QtCore import QRectF, Qt
from PyQt6.QtGui import QIcon, QTransform
from PyQt6.QtWidgets import (
    QGraphicsView,
    QHBoxLayout,
    QMainWindow,
    QPushButton,
    QVBoxLayout,
    QWidget,
)

from chi_editor.canvas import Canvas
from chi_editor.constants import ASSETS
from chi_editor.toolbar import CanvasToolBar

from .tasks.task import Task
from .bases.molecule.molecule import Molecule
from .bases.alpha_atom import AlphaAtom
from .playground import mol_from_graphs
class Editor(QMainWindow):
    # Hierarchy:
    #
    #   window:
    #       workspace:
    #           graphics_view:
    #               canvas:
    #       toolbar:

    # Contains everything except toolbar

    workspace: "QWidget"

    # QGraphicsView contains drawable space
    graphics_view: "QGraphicsView"

    # GraphicsScene where to draw all graphical objects
    canvas: "Canvas"

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

        # Add toolbar
        for index in range(3):
            self.toolbars[index] = CanvasToolBar(canvas=self.canvases[index], parent=self)
            self.addToolBar(Qt.ToolBarArea.LeftToolBarArea, self.toolbars[index])
        for toolbar in self.toolbars:
            toolbar.toggleViewAction().setChecked(False)
            toolbar.toggleViewAction().trigger()
            toolbar.toggleViewAction().trigger()
        self.toolbars[0].toggleViewAction().trigger()

        # Create dialogs (they won't show for now)
        self.choose_task_dialog = ChooseTaskDialog(editor=self)
        self.choose_task_dialog.setModal(True)

        self.result_dialog = TaskResultDialog(editor=self)
        self.result_dialog.setModal(True)

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

        solve_layout.addWidget(submit)
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

    def submitAnswer(self) -> None:
        items = self.canvases[EditorMode.SOLVE_MODE.value].items()

        if len(items) == 0:
            self.openResultDialog("No answer found")
            return


    def zoom_in(self) -> "None":
        # Get the current scale factor of the view
        current_scale = self.graphics_view.transform().m11()

        # Простите пожалуйста, Please forgive me, Entschuldigung bitte
        for item in items:
            if isinstance(item, AlphaAtom):
                molecule = item.molecule
                break

    def zoom_out(self) -> "None":
        # Get the current scale factor of the view
        current_scale = self.graphics_view.transform().m11()

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
