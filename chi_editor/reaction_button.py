from PyQt6.QtGui import QAction, QIcon


class ReactionButton(QAction):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setText("Reaction")
        self.setIcon(QIcon("..\\resources\\reaction.png"))
