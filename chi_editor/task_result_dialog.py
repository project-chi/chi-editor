from typing import TYPE_CHECKING

from PyQt6.QtWidgets import QDialog, QLabel, QSizePolicy, QVBoxLayout

if TYPE_CHECKING:
    from .editor import Editor


class TaskResultDialog(QDialog):
    # Main window
    editor: "Editor"

    # Text label
    _label: QLabel

    def __init__(self, *args, editor: "Editor", **kwargs):
        super().__init__(*args, **kwargs)
        self.editor = editor
        self._label = QLabel(self)

        layout = QVBoxLayout(self)
        layout.addWidget(self._label)
        self._label.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)

    def setText(self, text: str) -> None:
        self._label.setText(text)
