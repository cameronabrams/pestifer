# Author: ChatGPT 5
"""
A simple unified diff parser and applier.
"""

import re
from typing import List, Tuple

class PatchApplyError(Exception):
    pass

_hunk_re = re.compile(r'^@@ -(\d+)(?:,(\d+))? \+(\d+)(?:,(\d+))? @@')

def _parse_unified_hunks(patch_text: str) -> List[Tuple[int,int,int,int,List[str]]]:
    """
    Returns a list of hunks as tuples:
      (old_start, old_count, new_start, new_count, hunk_lines)
    where hunk_lines include the leading ' ', '+', '-' markers and keep newlines.
    """
    lines = patch_text.splitlines(keepends=True)
    hunks = []
    i = 0

    # Skip file header lines ('---', '+++') and any preamble
    while i < len(lines) and not lines[i].startswith('@@'):
        i += 1

    while i < len(lines):
        m = _hunk_re.match(lines[i])
        if not m:
            # allow interspersed non-hunk lines (e.g., additional headers)
            i += 1
            continue

        old_start = int(m.group(1))
        old_count = int(m.group(2) or '1')
        new_start = int(m.group(3))
        new_count = int(m.group(4) or '1')
        i += 1

        hunk_lines: List[str] = []
        # Collect lines until next '@@' or end. We can't rely only on counts
        # because context lines are mixed with additions/deletions.
        while i < len(lines) and not lines[i].startswith('@@'):
            # Stop when a new file header begins (rare inside a single-file patch)
            if lines[i].startswith(('--- ', '+++ ')) and hunk_lines:
                break
            hunk_lines.append(lines[i])
            i += 1

        hunks.append((old_start, old_count, new_start, new_count, hunk_lines))

    if not hunks:
        raise PatchApplyError("No unified diff hunks ('@@ ... @@') found. "
                              "Make sure this is a unified diff (use `diff -u`).")
    return hunks

def apply_unified_diff(content: str, patch_text: str, *, reverse: bool=False) -> str:
    """
    Apply a unified diff to `content` and return the patched string.

    - Set `reverse=True` to apply the patch in reverse (like `patch -R`).
    - Raises PatchApplyError if context does not match.
    """
    # Normalize to line-wise operations while preserving exact line endings
    src_lines = content.splitlines(keepends=True)
    out_lines: List[str] = []
    src_idx = 0  # 0-based index in src_lines

    hunks = _parse_unified_hunks(patch_text)

    for (old_start, old_count, new_start, new_count, hunk_lines) in hunks:
        # Unified diffs use 1-based line numbers for the OLD file
        target_idx = old_start - 1

        # Copy through all unchanged lines before this hunk
        if target_idx < src_idx:
            raise PatchApplyError(
                f"Hunk starting at old line {old_start} overlaps earlier edits."
            )
        out_lines.extend(src_lines[src_idx:target_idx])
        src_idx = target_idx

        # Apply hunk body
        for raw in hunk_lines:
            if raw.startswith('\\ No newline at end of file'):
                # meta line; ignore
                continue

            if reverse:
                # Swap meanings for reverse apply
                head = raw[0]
                if head == '+': head = '-'
                elif head == '-': head = '+'
                line = head + raw[1:]
            else:
                line = raw

            tag = line[:1]
            text = line[1:]

            if tag == ' ':
                # Context: must match input exactly
                if src_idx >= len(src_lines) or src_lines[src_idx] != text:
                    got = src_lines[src_idx] if src_idx < len(src_lines) else '<EOF>'
                    raise PatchApplyError(
                        "Context mismatch applying hunk.\n"
                        f"Expected: {text!r}\nGot     : {got!r}"
                    )
                out_lines.append(src_lines[src_idx])
                src_idx += 1

            elif tag == '-':
                # Deletion: input must match; do not copy to output
                if src_idx >= len(src_lines) or src_lines[src_idx] != text:
                    got = src_lines[src_idx] if src_idx < len(src_lines) else '<EOF>'
                    raise PatchApplyError(
                        "Deletion mismatch applying hunk.\n"
                        f"Expected to remove: {text!r}\nBut found          : {got!r}"
                    )
                src_idx += 1

            elif tag == '+':
                # Addition: copy to output, do not consume input
                out_lines.append(text)

            elif tag in ('@', '-', '+'):
                # Should have been handled; included for completeness
                raise PatchApplyError(f"Unexpected hunk line: {line!r}")

            else:
                # Unknown prefix (could be stray headers); be strict
                raise PatchApplyError(f"Invalid hunk line (no prefix): {line!r}")

    # Copy any remaining source lines
    out_lines.extend(src_lines[src_idx:])
    return ''.join(out_lines)
