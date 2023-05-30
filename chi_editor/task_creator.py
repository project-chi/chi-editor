from PyQt6.QtWidgets import QDialog, QLineEdit, QDialogButtonBox, QFormLayout, QComboBox
from PyQt6.QtCore import Qt, pyqtSignal

from chi_editor.api.server import Server, default_url
from chi_editor.api.task import Kind
from chi_editor.canvas import Canvas
from chi_editor.chains.chain import Chain
from chi_editor.reactions.reaction import Reaction
from chi_editor.utils.json_utils import create_task
from chi_editor.editor_mode import EditorMode


class InputDialog(QDialog):
    canvas_molecule: Canvas
    canvas_reaction: Canvas
    canvas_chain: Canvas
    active_canvas: Canvas
    canvas_changed: pyqtSignal = pyqtSignal(Canvas, int)

    def __init__(self, canvas_molecule: Canvas, canvas_reaction: Canvas, canvas_chain: Canvas, parent=None):
        super().__init__(parent)

        self.canvas_molecule = canvas_molecule
        self.canvas_reaction = canvas_reaction
        self.canvas_chain = canvas_chain
        self.active_canvas = canvas_molecule
        canvas_molecule.canvas_in_focus.connect(lambda: self.changeCanvas(1))
        canvas_chain.canvas_in_focus.connect(lambda: self.changeCanvas(0))
        canvas_reaction.canvas_in_focus.connect(lambda: self.changeCanvas(2))
        self.name = QLineEdit(self)
        self.formulation = QLineEdit(self)
        self.correct_answer = QLineEdit(self)
        self.type = QComboBox(self)
        button_box = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok | \
                                      QDialogButtonBox.StandardButton.Help | \
                                      QDialogButtonBox.StandardButton.Cancel | \
                                      QDialogButtonBox.StandardButton.Reset, self)
        button_box.button(QDialogButtonBox.StandardButton.Help).setText("Parse")
        button_box.button(QDialogButtonBox.StandardButton.Reset).setText("Create Local")
        # Fill type combo box
        for kind in Kind:
            self.type.addItem(kind.name, userData=kind)
        self.type.setCurrentIndex(1)
        self.type.currentIndexChanged.connect(self.changeCanvas)

        layout = QFormLayout(self)
        layout.addRow("Choose type of the task", self.type)
        layout.addRow("Name of the task", self.name)
        layout.addRow("Formulate the problem", self.formulation)
        layout.addRow("Input the correct answer as SMILES", self.correct_answer)
        layout.addWidget(button_box)

        button_box.accepted.connect(self.sendTask)
        button_box.rejected.connect(self.clearAll)
        button_box.helpRequested.connect(self.parseInput)

        button_box.button(QDialogButtonBox.StandardButton.Reset).clicked.connect(self.createTask)

    def changeCanvas(self, index: int):
        match index:
            case 0:
                self.active_canvas = self.canvas_chain
            case 1:
                self.active_canvas = self.canvas_molecule
            case 2:
                self.active_canvas = self.canvas_reaction
        self.canvas_changed.emit(self.active_canvas, EditorMode.CREATE_MODE.value)

    def getInputs(self):
        current_type = self.type.currentData(Qt.ItemDataRole.UserRole)
        return self.name.text(), current_type, self.formulation.text(), self.correct_answer.text()

    def clearAll(self):
        self.name.clear()
        self.formulation.clear()
        self.correct_answer.clear()

    def sendTask(self):
        res_name, res_type, res_formulation, res_correct = self.getInputs()
        server = Server(default_url)
        server.create_task(res_name, res_type, res_formulation, res_correct)
        self.clearAll()

    def parseInput(self):
        answer_type: type
        if self.active_canvas == self.canvas_molecule:
            self.correct_answer.setText(self.active_canvas.findMolecule())
            return
        elif self.active_canvas == self.canvas_reaction:
            answer_type = Reaction
        else:
            answer_type = Chain
        self.correct_answer.setText(self.active_canvas.findElement(answer_type))

    def createTask(self):
        res_name, res_type, res_formulation, res_correct = self.getInputs()
        create_task(res_name, res_type.value, res_formulation, res_correct, "")
        self.clearAll()
