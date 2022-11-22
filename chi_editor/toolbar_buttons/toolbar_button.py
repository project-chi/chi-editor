from PyQt6.QtGui import QAction


class ToolBarButton(QAction):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setCheckable(True)
