"""Where a design job and its results live between the submit and the reader.

A run takes minutes and the browser tab that started it is long gone by the end, so the job is a
record on disk rather than something held in memory: a row in SQLite for the status and the
parameters, and the result tables beside it under the data directory. That is what lets the results
be opened later, from a link, by whoever has it.
"""

import json
import os
import sqlite3
import uuid
from datetime import datetime, timezone
from pathlib import Path

QUEUED, RUNNING, DONE, FAILED = "queued", "running", "done", "failed"

# The tables a finished job leaves behind, in the order the results page shows them.
RESULT_FILES = ["designed_asos.csv", "safety_detail.csv", "off_targets.csv"]


def jobs_dir() -> Path:
    root = Path(os.environ.get("TAUSO_DATA_DIR", "/home/mambauser/.tauso_data")) / "jobs"
    root.mkdir(parents=True, exist_ok=True)
    return root


def job_dir(job_id: str) -> Path:
    return jobs_dir() / job_id


def _connect():
    # Several processes touch this: the UI reads while a worker writes, so WAL rather than the
    # default rollback journal, which would lock readers out for the length of a write.
    connection = sqlite3.connect(jobs_dir() / "jobs.db", timeout=30)
    connection.row_factory = sqlite3.Row
    connection.execute("PRAGMA journal_mode=WAL")
    connection.execute(
        """CREATE TABLE IF NOT EXISTS jobs (
               id TEXT PRIMARY KEY,
               status TEXT NOT NULL,
               created_at TEXT NOT NULL,
               finished_at TEXT,
               target TEXT,
               source_info TEXT,
               email TEXT,
               parameters TEXT,
               error TEXT
           )"""
    )
    return connection


def _now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def create(target: str, source_info: str, email: str, parameters: dict) -> str:
    """Record a job as queued and return the id its results will be found under."""
    job_id = uuid.uuid4().hex[:12]
    job_dir(job_id).mkdir(parents=True, exist_ok=True)
    with _connect() as connection:
        connection.execute(
            "INSERT INTO jobs (id, status, created_at, target, source_info, email, parameters)"
            " VALUES (?, ?, ?, ?, ?, ?, ?)",
            (job_id, QUEUED, _now(), target, source_info, email, json.dumps(parameters)),
        )
    return job_id


def mark(job_id: str, status: str, error: str = None) -> None:
    finished = _now() if status in (DONE, FAILED) else None
    with _connect() as connection:
        connection.execute(
            "UPDATE jobs SET status = ?, finished_at = COALESCE(?, finished_at), error = ? WHERE id = ?",
            (status, finished, error, job_id),
        )


def get(job_id: str) -> dict:
    """The job with this id, or None. `parameters` comes back as a dict."""
    with _connect() as connection:
        row = connection.execute("SELECT * FROM jobs WHERE id = ?", (job_id,)).fetchone()
    if row is None:
        return None
    job = dict(row)
    job["parameters"] = json.loads(job["parameters"] or "{}")
    return job


def recent(limit: int = 20) -> list:
    with _connect() as connection:
        rows = connection.execute(
            "SELECT id, status, created_at, target FROM jobs ORDER BY created_at DESC LIMIT ?", (limit,)
        ).fetchall()
    return [dict(row) for row in rows]


def save_results(job_id: str, tables: dict) -> None:
    """Write the result tables for a job, named as RESULT_FILES expects to find them."""
    directory = job_dir(job_id)
    directory.mkdir(parents=True, exist_ok=True)
    for name, frame in tables.items():
        frame.to_csv(directory / name, index=False)


def results_path(job_id: str, name: str) -> Path:
    return job_dir(job_id) / name


def has_results(job_id: str) -> bool:
    return all(results_path(job_id, name).exists() for name in RESULT_FILES)


def public_url(job_id: str) -> str:
    """The address to hand someone so they can open these results."""
    base = os.environ.get("PUBLIC_BASE_URL", "http://localhost:8501").rstrip("/")
    return f"{base}/?job={job_id}"
