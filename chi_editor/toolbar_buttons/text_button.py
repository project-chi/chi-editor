from .toolbar_button import ToolBarButton


class TextButton(ToolBarButton):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setText("Text")
