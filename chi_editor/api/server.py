from uuid import UUID

import requests

from chi_editor.api.task import Kind, Task

default_url: str = "https://project-chi.kapkekes.site/api"

class Server:
    _server_url: str

    def __init__(self, server_url: str) -> None:
        swagger_response = requests.get(f"{server_url}/documentation")
        if not swagger_response.ok:
            raise ValueError("can't find server at this URL")

        self._server_url = server_url

    def get_tasks(self) -> list[UUID]:
        request_url = f"{self._server_url}/tasks"
        response = requests.get(request_url)
        return [UUID(identifier) for identifier in response.json()]

    def create_task(
        self,
        name: str,
        kind: Kind,
        problem: str,
        solution: str,
        initial: str = "",
    ) -> Task:
        request_url = f"{self._server_url}/tasks"
        payload = {
            "name": name,
            "kind": kind.value,
            "problem": problem,
            "solution": solution,
            "initial": initial,
        }
        response = requests.post(request_url, json=payload)
        if not response.ok:
            raise ValueError(response.json())

        return Task.parse_obj(response.json())

    def get_task(self, uuid: str | UUID) -> Task:
        if isinstance(uuid, str):
            uuid = UUID(uuid)

        request_url = f"{self._server_url}/tasks/{uuid}"
        response = requests.get(request_url)

        if not response.ok:
            raise ValueError(f"there is no post with {uuid = }")

        return Task.parse_obj(response.json())

    def update_task(
        self,
        uuid: str | UUID,
        name: str | None = None,
        kind: Kind | None = None,
        problem: str | None = None,
        solution: str | None = None,
        initial: str | None = None,
    ) -> Task:
        if isinstance(uuid, str):
            uuid = UUID(uuid)

        request_url = f"{self._server_url}/tasks/{uuid}"
        payload = {
            "name": name,
            "kind": kind.value if kind is not None else None,
            "problem": problem,
            "solution": solution,
            "initial": initial,
        }
        payload = {
            key: val for key, val in payload.items()
            if val is not None
        }
        response = requests.patch(request_url, json=payload)
        if not response.ok:
            raise ValueError(f"there is no post with {uuid = }")

        return Task.parse_obj(response.json())

    def delete_task(self, uuid: str | UUID) -> Task:
        if isinstance(uuid, str):
            uuid = UUID(uuid)

        request_url = f"{self._server_url}/tasks/{uuid}"
        response = requests.delete(request_url)

        if not response.ok:
            raise ValueError(f"there is no post with {uuid = }")

        return Task.parse_obj(response.json())
