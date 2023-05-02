from PyQt6.QtWidgets import (
    QDialog,
    QLineEdit,
    QDialogButtonBox,
    QFormLayout,
    QComboBox,
)
from PyQt6.QtCore import Qt

from chi_editor.api.server import Server, default_url
from chi_editor.api.task import Kind
from chi_editor.canvas import Canvas
from chi_editor.utils.json_utils import create_task


class InputDialog(QDialog):
    canvas: Canvas

    def __init__(self, canvas: Canvas, parent=None):
        super().__init__(parent)

        self.canvas = canvas
        self.name = QLineEdit(self)
        self.formulation = QLineEdit(self)
        self.correct_answer = QLineEdit(self)
        self.type = QComboBox(self)
        button_box = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok
            | QDialogButtonBox.StandardButton.Help
            | QDialogButtonBox.StandardButton.Cancel
            | QDialogButtonBox.StandardButton.Reset,
        )
        button_box.button(QDialogButtonBox.StandardButton.Help).setText("Parse")
        button_box.button(QDialogButtonBox.StandardButton.Reset).setText("Create Local")
        # Fill type combo box
        for kind in Kind:
            self.type.addItem(kind.name, userData=kind)

        layout = QFormLayout(self)
        layout.addRow("Choose type of the task", self.type)
        layout.addRow("Name of the task", self.name)
        layout.addRow("Formulate the problem", self.formulation)
        layout.addRow("Input the correct answer as SMILES", self.correct_answer)
        layout.addWidget(button_box)

        button_box.accepted.connect(self.sendTask)
        button_box.rejected.connect(self.clearAll)
        button_box.helpRequested.connect(self.parseInput)

        button_box.button(QDialogButtonBox.StandardButton.Reset).clicked.connect(
            self.createTask
        )

    def getInputs(self):
        current_type = self.type.currentData(Qt.ItemDataRole.UserRole)
        return (
            self.name.text(),
            current_type,
            self.formulation.text(),
            self.correct_answer.text(),
        )

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
        self.correct_answer.setText(self.canvas.findMolecule())

    def createTask(self):
        res_name, res_type, res_formulation, res_correct = self.getInputs()
        create_task(res_name, res_type.value, res_formulation, res_correct, "")
        self.clearAll()
