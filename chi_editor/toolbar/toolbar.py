from chi_editor.toolbar.general_toolbar import GeneralToolBar
from chi_editor.toolbar.tools.answer_field_tool import AnswerFieldTool


class CanvasToolBar(GeneralToolBar):
    def __init__(self, *args, **kwargs) -> "None":
        super().__init__(*args, **kwargs)
        tool = AnswerFieldTool(self._canvas)
        self.addAction(tool)
        self.action_group.addAction(tool)
        self.widgetForAction(tool).setStyleSheet("padding: 5px")
