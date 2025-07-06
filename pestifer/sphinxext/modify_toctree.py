"""  
Modifies toctree directives in documentation when ``pestifer modify-package``` is used
to add or remove examples.  Written substantially by ChatGPT, with some modifications.
"""

import os

def detect_common_prefix(entries):
    split_entries = [e.split(os.sep) for e in entries if os.sep in e]
    if not split_entries:
        return None
    first_dir = split_entries[0][0]
    if all(parts[0] == first_dir for parts in split_entries):
        return first_dir + os.sep
    return None

def read_rst_file(filepath):
    with open(filepath, "r", encoding="utf-8") as f:
        return f.readlines()

def write_rst_file(filepath, lines):
    with open(filepath, "w", encoding="utf-8") as f:
        f.writelines(lines)

def find_toctree_block(lines):
    for i, line in enumerate(lines):
        if line.strip().startswith(".. toctree::"):
            start = i
            break
    else:
        raise ValueError("No toctree block found.")

    # Collect options
    i = start + 1
    while i < len(lines) and lines[i].lstrip().startswith(":"):
        i += 1
    entry_start = i

    # Collect entries (non-empty lines starting with indentation)
    entry_end = entry_start
    while entry_end < len(lines) and (lines[entry_end].startswith("  ") or lines[entry_end].strip() == ""):
        entry_end += 1

    return start, entry_start, entry_end

def parse_toctree_entries(lines, entry_start, entry_end):
    return [line.strip() for line in lines[entry_start:entry_end] if line.strip()]

def reconstruct_toctree_block(lines, start, entry_start, entry_end, new_entries):
    # Keep the toctree line and any options
    header = lines[start:entry_start]

    # Ensure there's exactly one blank line between options and entries
    if not header or not header[-1].strip() == "":
        header.append("\n")

    # Format entries with two-space indent
    formatted_entries = [f"  {entry}\n" for entry in new_entries]

    return lines[:start] + header + formatted_entries + lines[entry_end:]

def modify_entries(entries, action, target=None, new_entry=None):
    prefix = detect_common_prefix(entries)
    entries_set = set(entries)

    def apply_prefix(e):
        if os.sep in e or not prefix:
            return e
        return prefix + e

    if action == "add":
        entry_to_add = apply_prefix(new_entry)
        if entry_to_add not in entries_set:
            entries.append(entry_to_add)

    elif action == "delete":
        entries = [e for e in entries if e != target]

    elif action == "insert":
        if target in entries:
            idx = entries.index(target)
            entry_to_insert = apply_prefix(new_entry)
            if entry_to_insert not in entries:
                entries.insert(idx + 1, entry_to_insert)

    return entries


def modify_toctree(filepath, action, target=None, new_entry=None):
    lines = read_rst_file(filepath)
    start, entry_start, entry_end = find_toctree_block(lines)
    entries = parse_toctree_entries(lines, entry_start, entry_end)
    updated_entries = modify_entries(entries, action, target, new_entry)
    new_lines = reconstruct_toctree_block(lines, start, entry_start, entry_end, updated_entries)
    write_rst_file(filepath, new_lines)

# Examples:
# update_rst("examples.rst", action="add", new_entry="examples/new_example")
# update_rst("examples.rst", action="delete", target="examples/8fad")
# update_rst("examples.rst", action="insert", target="examples/env", new_entry="examples/new_between")

