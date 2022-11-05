from PyQt6.QtWidgets import QMenu


class HelpMenu(QMenu):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setTitle("Help")
