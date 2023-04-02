from PyQt6.QtWidgets import QDialog, QLineEdit, QDialogButtonBox, QFormLayout


class InputDialog(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)

        self.name = QLineEdit(self)
        self.formulation = QLineEdit(self)
        self.correct_answer = QLineEdit(self)
        self.type = QLineEdit(self)
        button_box = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel, self)

        layout = QFormLayout(self)
        layout.addRow("Name of the task", self.name)
        layout.addRow("Formulate the problem", self.formulation)
        layout.addRow("Input the correct answer as SMILES", self.correct_answer)
        layout.addRow("What type is this problem of?", self.type)
        layout.addWidget(button_box)

        button_box.accepted.connect(self.accept)
        button_box.rejected.connect(self.reject)

    def getInputs(self):
        return self.first.text(), self.second.text()
