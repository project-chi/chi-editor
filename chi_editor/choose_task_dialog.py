from PyQt6.QtWidgets import QDialog


class ChooseTaskDialog(QDialog):
    def __init__(self, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs)
        self.setWindowTitle("Choose a task")
