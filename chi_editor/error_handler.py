from PyQt6.QtWidgets import QMessageBox


def error_handler(error_message: str) -> None:
    error_box = QMessageBox()
    error_box.setIcon(QMessageBox.Icon.Critical)
    error_box.setWindowTitle("Oopsie Doopsie")
    error_box.setText(error_message)
    error_box.setStandardButtons(QMessageBox.StandardButton.Ok)
    error_box.exec()
