from rdkit import Chem

from chi_editor.tasks.task_type import TaskType


class Task:
    # Task's title
    _title: str

    # Answer represented as one SMILES string
    _correct_answer: str

    # Task's type
    _type: TaskType

    def __init__(self, title: str, formulation: str, correct_answer: str, task_type: TaskType):
        self._title = title
        self._formulation = formulation
        self._correct_answer = Chem.CanonSmiles(correct_answer)
        self._type = task_type

    def checkAnswer(self, user_answer: str) -> bool:
        canon_user_answer = Chem.CanonSmiles(user_answer)
        return self._correct_answer == canon_user_answer

    def title(self) -> str:
        return self._title

    def correct_answer(self) -> str:
        return self._correct_answer

    def task_type(self) -> TaskType:
        return self._type

    def formulation(self) -> str:
        return self._formulation
    