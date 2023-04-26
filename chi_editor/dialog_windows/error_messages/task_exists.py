from PyQt6.QtWidgets import QMessageBox


class TaskExistsMessageBox(QMessageBox):
    def __init__(self, *args, **kwargs) -> None:
        super().__init__()
        self.setText("Task with this name already exists")
        self.setIcon(QMessageBox.Icon.Critical)
