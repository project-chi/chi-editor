from PyQt6.QtWidgets import QMenu


class ViewMenu(QMenu):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.setTitle("View")
