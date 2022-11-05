from PyQt6.QtGui import QAction


class ReactionButton(QAction):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setText("Reaction")
