from PyQt6.QtWidgets import QDialog, QLabel, QVBoxLayout, QSizePolicy


class TaskCreationDialog(QDialog):
    # Text label
    _label: QLabel

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._label = QLabel(self)
        self.setWindowTitle("Result")

        layout = QVBoxLayout(self)
        layout.addWidget(self._label)
        self._label.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)

    def setText(self, text: str) -> None:
        self._label.setText(text)
