from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import QGraphicsSceneMouseEvent, QInputDialog, QWidget
from rdkit import Chem

from .structure import put_molecule
from ...bases.tool import Tool


class Smiles(Tool):
    def mouse_press_event(self, event: QGraphicsSceneMouseEvent) -> None:
        if event.button() == Qt.MouseButton.LeftButton:
            dialog = SmilesDialog()
            self.canvas.addWidget(dialog)
            if dialog.smiles != "":
                molecule: Chem.Mol = Chem.MolFromSmiles(dialog.smiles)
                Chem.Kekulize(molecule)
                put_molecule(self.canvas, molecule, event.scenePos())
            self.canvas.removeItem(dialog.graphicsProxyWidget())

    @property
    def asset(self) -> str:
        return "structure"


class SmilesDialog(QWidget):
    smiles: str

    def __init__(self):
        super().__init__()
        self.initUI()

    def initUI(self):
        self.showDialog()

    def showDialog(self):
        text, ok = QInputDialog.getText(self, 'input dialog', 'Is this ok?')
        if ok:
            self.smiles = str(text)
        else:
            self.smiles = ""
