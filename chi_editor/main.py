import sys

from PyQt6.QtWidgets import QApplication, QMessageBox

from chi_editor.editor import Editor


def error_handler(error_message: str) -> None:
    error_box = QMessageBox()
    error_box.setIcon(QMessageBox.Icon.Critical)
    error_box.setWindowTitle("Oopsie Doopsie")
    error_box.setText(error_message)
    error_box.setStandardButtons(QMessageBox.StandardButton.Ok)
    error_box.exec()
    pass


def main():
    application = QApplication(sys.argv)
    window = Editor()
    window.show()
    sys.exit(application.exec())
