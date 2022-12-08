from PyQt6.QtWidgets import QMenu


class HelpMenu(QMenu):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.setTitle("Help")
