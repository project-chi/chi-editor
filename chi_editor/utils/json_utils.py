import json


def create_task(name: str, kind: int, problem: str, initial: str, solution: str):
    dictionary = {
        "name": name,
        "kind": kind,
        "problem": problem,
        "initial": initial,
        "solution": solution,
        "identifier": ""
    }

    json_object = json.dumps(dictionary)

    with open(name + ".json", "w") as outfile:
        outfile.write(json_object)
