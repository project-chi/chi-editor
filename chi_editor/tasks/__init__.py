from chi_editor.tasks.task import Task
from chi_editor.tasks.task_type import TaskType

alkanes_list: tuple[list[Task], TaskType] = (
    [
        Task("Methane", "C", TaskType.ALKANES),
        Task("Ethane", "CC", TaskType.ALKANES)
    ],
    TaskType.ALKANES
)

alcohols_list: tuple[list[Task], TaskType] = (
    [
        Task("Ethanol", "CCO", TaskType.ALCOHOLS),
        Task("1-Butanol", "CCCCO", TaskType.ALCOHOLS),
    ],
    TaskType.ALCOHOLS
)

aromatic_list: tuple[list[Task], TaskType] = (
    [
        Task("Benzene", "c1ccccc1", TaskType.AROMATIC),
        Task("Toluene", "Cc1ccccc1", TaskType.AROMATIC),
    ],
    TaskType.AROMATIC
)

tasks_list: list[tuple[list[Task], TaskType]] = [
    alkanes_list,
    alcohols_list,
    aromatic_list,
]
