from chi_editor.toolbar.general_toolbar import GeneralToolBar
from chi_editor.toolbar.tools.extra_tools import main_canvas_tools


class CanvasToolBar(GeneralToolBar):

    def __init__(self, *args, **kwargs) -> "None":
        super().__init__(*args, **kwargs)

        for Tool in main_canvas_tools:
            tool = Tool(self._canvas)
            self.addAction(tool)
            self.action_group.addAction(tool)
            self.widgetForAction(tool).setStyleSheet("padding: 5px")
