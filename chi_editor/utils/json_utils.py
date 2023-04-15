import json
from typing import ClassVar
from pathlib import Path

from chi_editor.constants import RESOURCES

# Default folder for local task files
default_dir: ClassVar[Path] = RESOURCES / "local_tasks"


def create_task(name: str, kind: int, problem: str, solution: str, initial: str):
    dictionary = {
        "name": name,
        "kind": kind,
        "problem": problem,
        "initial": initial,
        "solution": solution,
        "identifier": "00000000-0000-0000-0000-000000000000"
    }

    json_object = json.dumps(dictionary, indent=2)

    with open(str(default_dir.absolute() / (name + ".json")), "w") as outfile:
        outfile.write(json_object)
