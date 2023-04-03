from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import QGraphicsSceneMouseEvent, QInputDialog, QWidget
from rdkit import Chem

from ...bases.tool import Tool
from .structure import put_molecule


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
    def picture(self) -> str:
        return "smiles"


class SmilesDialog(QWidget):
    smiles: str

    def __init__(self):
        super().__init__()
        self.initUI()

    def initUI(self):
        self.showDialog()

    def showDialog(self):
        text, ok = QInputDialog.getText(
            self, "input dialog", "Put in your SMILES formula"
        )
        if ok:
            self.smiles = str(text)
        else:
            self.smiles = ""
