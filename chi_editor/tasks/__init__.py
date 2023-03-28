from task import Task
from task_type import TaskType

tasks_list: list[Task] = [
    Task("Methane", "C", TaskType.ALKANES),
    Task("Ethane", "CC", TaskType.ALKANES),
    Task("Ethanol", "CCO", TaskType.ALCOHOLS),
    Task("1-Butanol", "CCCCO", TaskType.ALCOHOLS),
    Task("Benzene", "c1ccccc1", TaskType.AROMATIC),
    Task("Toluene", "Cc1ccccc1", TaskType.AROMATIC)
]
