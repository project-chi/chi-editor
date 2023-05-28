from chi_editor.bases.tool import Tool
from chi_editor.toolbar.tools.answer_field_tool import AnswerFieldTool
from chi_editor.toolbar.tools.add_reaction import AddReaction
from chi_editor.toolbar.tools.add_chain import AddChain

main_canvas_tools: tuple[type[Tool], ...] = (
    AnswerFieldTool,
    AddReaction,
    AddChain
)
