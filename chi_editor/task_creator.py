from PyQt6.QtWidgets import QDialog, QLineEdit, QDialogButtonBox, QFormLayout, QComboBox

from chi_editor.api.server import Server
from chi_editor.api.task import Kind


class InputDialog(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)

        self.name = QLineEdit(self)
        self.formulation = QLineEdit(self)
        self.correct_answer = QLineEdit(self)
        self.type = QComboBox(self)
        button_box = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel, self)

        # Fill type combo box
        for kind in Kind:
            self.type.addItem(kind.name)

        layout = QFormLayout(self)
        layout.addRow("Choose type of the task", self.type)
        layout.addRow("Name of the task", self.name)
        layout.addRow("Formulate the problem", self.formulation)
        layout.addRow("Input the correct answer as SMILES", self.correct_answer)
        layout.addWidget(button_box)

        button_box.accepted.connect(self.sendTask)
        button_box.rejected.connect(self.reject)

    def getInputs(self):
        return self.name.text(), self.formulation.text(), self.correct_answer.text()

    def sendTask(self):
        res_name, res_formulation, res_correct = self.getInputs()
        server = Server("http://kapkekes.site:8000")
        server.create_task(res_name, Kind.Molecule, res_formulation, res_correct)
