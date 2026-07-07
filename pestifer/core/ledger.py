# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
The **modification ledger**: an append-only, git-tracked record of the changes that
``pestifer modify-package`` applies to the package (custom residues, PDB-repository entries,
examples, force-field regenerations).  It complements the resource *provenance* tags shown by
``show-resources`` -- provenance says a file is ``custom``; the ledger says *who added it,
when, via which operation, and on what branch*, and lets a modification be reverted.

The ledger lives at ``pestifer/resources/modifications.jsonl`` -- one JSON object per line, so
each contribution appends a single line (a clean, reviewable diff and minimal merge conflicts).
"""
import json

from datetime import datetime, timezone
from pathlib import Path

LEDGER_FILENAME = 'modifications.jsonl'
SCHEMA_VERSION = 1


def ledger_file(resources_path) -> Path:
    """Path to the ledger file under ``resources_path``."""
    return Path(resources_path) / LEDGER_FILENAME


def read(resources_path) -> list:
    """Return all ledger entries (oldest first); an empty list if the ledger does not exist."""
    path = ledger_file(resources_path)
    if not path.is_file():
        return []
    entries = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line:
                entries.append(json.loads(line))
    return entries


def get(resources_path, entry_id: int):
    """Return the entry with ``id == entry_id``, or ``None``."""
    for e in read(resources_path):
        if e.get('id') == entry_id:
            return e
    return None


def _next_id(entries) -> int:
    return max((e.get('id', 0) for e in entries), default=0) + 1


def _now() -> str:
    """Current UTC time as an ISO-8601 string (seconds precision)."""
    return datetime.now(timezone.utc).isoformat(timespec='seconds')


def append(resources_path, *, category, verb, summary, files=None, author='',
           branch=None, extra=None, timestamp=None) -> dict:
    """
    Append a new ledger entry and return it (with its assigned ``id``).  ``files`` should be
    repo-relative paths.  ``timestamp`` may be supplied for reproducible tests; otherwise the
    current UTC time is used.
    """
    entries = read(resources_path)
    entry = {
        'schema': SCHEMA_VERSION,
        'id': _next_id(entries),
        'timestamp': timestamp or _now(),
        'category': category,
        'verb': verb,
        'summary': summary,
        'files': list(files or []),
        'author': author,
        'branch': branch,
        'reverted': False,
    }
    if extra:
        entry.update(extra)
    with open(ledger_file(resources_path), 'a') as f:
        f.write(json.dumps(entry) + '\n')
    return entry


def write_all(resources_path, entries):
    """Rewrite the whole ledger (used to update an entry's status, e.g. mark it reverted)."""
    with open(ledger_file(resources_path), 'w') as f:
        for e in entries:
            f.write(json.dumps(e) + '\n')


def format_entries(entries, limit: int = None, category: str = None) -> str:
    """Render ledger entries as a compact, human-readable listing (newest last)."""
    rows = [e for e in entries if category is None or e.get('category') == category]
    if limit:
        rows = rows[-limit:]
    if not rows:
        return '(no ledger entries)'
    lines = []
    for e in rows:
        date = str(e.get('timestamp', ''))[:10]
        tag = f"#{e.get('id')}"
        status = ''
        if e.get('reverted'):
            by = e.get('reverted_by')
            status = f"  [reverted{f' by #{by}' if by else ''}]"
        elif e.get('reverts') is not None:
            status = f"  [reverts #{e.get('reverts')}]"
        lines.append(f"{tag:>5}  {date}  {e.get('category')}/{e.get('verb')}"
                     f"  {e.get('summary', '')}{status}")
    return '\n'.join(lines)
