"""  
Modifies toctree directives in documentation when ``pestifer modify-package`` is used
to add or remove examples.  Written substantially by ChatGPT, with some modifications.
"""

import os

def detect_common_prefix(entries):
    """
    Detects the common prefix in a list of entries.
    If all entries share a common prefix, returns that prefix with a trailing separator.
    If no common prefix exists, returns None.
    """
    split_entries = [e.split(os.sep) for e in entries if os.sep in e]
    if not split_entries:
        return None
    first_dir = split_entries[0][0]
    if all(parts[0] == first_dir for parts in split_entries):
        return first_dir + os.sep
    return None

def read_rst_file(filepath):
    """
    Reads a reStructuredText (RST) file and returns its lines.
    
    Parameters
    ----------
    filepath : str
        The path to the RST file to read.
        
    Returns
    -------
    list of str
        A list of lines from the RST file.
    """
    with open(filepath, "r", encoding="utf-8") as f:
        return f.readlines()

def write_rst_file(filepath, lines):
    """
    Writes a list of lines to a reStructuredText (RST) file.
    
    Parameters
    ----------
    filepath : str
        The path to the RST file to write.
    lines : list of str
        The lines to write to the RST file.
    """
    with open(filepath, "w", encoding="utf-8") as f:
        f.writelines(lines)

def find_toctree_block(lines):
    """
    Finds the start and end indices of the toctree block in the RST file lines.
    The toctree block starts with ".. toctree::" and may have options before the entries.

    Parameters
    ----------
    lines : list of str
        The lines of the RST file.

    Returns
    -------
    tuple
        A tuple containing the start index of the toctree block, the start index of the
        entries, and the end index of the entries.

    Raises
    ------
    ValueError
        If no toctree block is found in the lines.
    """
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
    """
    Parses the entries from the toctree block in the RST file lines.
    
    Parameters
    ----------
    lines : list of str
        The lines of the RST file.
    entry_start : int
        The start index of the entries in the toctree block.
    entry_end : int
        The end index of the entries in the toctree block.
  
    Returns
    -------
    list of str
        A list of entries from the toctree block, stripped of leading and trailing whitespace.
    """
    return [line.strip() for line in lines[entry_start:entry_end] if line.strip()]

def reconstruct_toctree_block(lines, start, entry_start, entry_end, new_entries):
    """
    Reconstructs the toctree block in the RST file lines with the updated entries
    
    Parameters
    ----------
    lines : list of str
        The original lines of the RST file.
    start : int
        The start index of the toctree block.
    entry_start : int
        The start index of the entries in the toctree block.
    entry_end : int
        The end index of the entries in the toctree block.
    new_entries : list of str
        The new entries to be included in the toctree block. 
        
    Returns
    -------
    list of str
        The lines of the RST file with the updated toctree block."""
    # Keep the toctree line and any options
    header = lines[start:entry_start]

    # Ensure there's exactly one blank line between options and entries
    if not header or not header[-1].strip() == "":
        header.append("\n")

    # Format entries with two-space indent
    formatted_entries = [f"  {entry}\n" for entry in new_entries]

    return lines[:start] + header + formatted_entries + lines[entry_end:]

def modify_entries(entries, action, target=None, new_entry=None, index=None, common_prefix=None):
    """
    Modifies the list of entries based on the specified action.
    
    Parameters
    ----------
    entries : list of str
        The current list of entries in the toctree.
    action : str
        The action to perform: "delete", "update", or "insert".
    target : str, optional
        The entry to delete or update if action is "delete" or "update", respectively. Required for "delete" and "update" actions.
    new_entry : str, optional
        The new entry to insert. Required for "insert" actions.   
    index : int, optional
        The 1-based index at which to insert the new entry. Required for "insert" and "update" actions.

    Returns
    -------
    list of str
        The modified list of entries in the toctree.
    """

    if action not in ["delete", "insert", "update"]:
        raise ValueError(f"Invalid action: {action}. Must be 'delete', 'insert', or 'update'.")
    prefix = detect_common_prefix(entries)
    if not prefix and common_prefix is not None:
        prefix = common_prefix+os.sep
    def apply_prefix(e):
        if os.sep in e or not prefix:
            return e
        return prefix + e

    if action == "delete":
        entries = [e for e in entries if e != apply_prefix(target)]

    elif action == "insert":
        entry_to_insert = apply_prefix(new_entry)
        if entry_to_insert not in entries:
            entries.insert(index-1, entry_to_insert)
    
    elif action== "update":
        if apply_prefix(target) in entries:
            entries[entries.index(apply_prefix(target))] = apply_prefix(new_entry)
        else:
            raise ValueError(f"Target entry {target} not found for update.")

    return entries

def modify_toctree(filepath, action, target=None, new_entry=None, index=None, common_prefix=None):
    """
    Modifies the toctree block in a reStructuredText (RST) file based on the specified action.
    
    Parameters
    ----------
    filepath : str
        The path to the RST file to modify.
    action : str
        The action to perform: "delete", "update", or "insert".
    target : str, optional
        The entry to modify (delete, update, or insert before/after). Required for "delete", "insert", and "update" actions.
    new_entry : str, optional
        The new entry to add or insert. Required for "insert" and "update" actions.
    index : int, optional
        The 1-based index at which to insert the new entry or update an existing entry. Required for "insert" and "update" actions.

    """
    lines = read_rst_file(filepath)
    start, entry_start, entry_end = find_toctree_block(lines)
    entries = parse_toctree_entries(lines, entry_start, entry_end)
    updated_entries = modify_entries(entries, action, target, new_entry, index, common_prefix=common_prefix)
    new_lines = reconstruct_toctree_block(lines, start, entry_start, entry_end, updated_entries)
    write_rst_file(filepath, new_lines)

def get_num_entries_in_toctree(filepath):
    """
    Retrieves the number of entries in the toctree block of an RST file.
    
    Parameters
    ----------
    filepath : str
        The path to the RST file.
        
    Returns
    -------
    int
        The number of entries in the toctree block.
        
    Raises
    ------
    ValueError
        If no toctree block is found in the file.
    """
    if not os.path.isfile(filepath):
        return 0
    lines = read_rst_file(filepath)
    start, entry_start, entry_end = find_toctree_block(lines)
    entries = parse_toctree_entries(lines, entry_start, entry_end)
    return len(entries)

def get_name_from_toctree(filepath, index):
    """
    Retrieves the name of the entry at the specified index from the toctree in an RST file.
    
    Parameters
    ----------
    filepath : str
        The path to the RST file.
    index : int
        The 1-based index of the entry to retrieve.
        
    Returns
    -------
    str
        The name of the entry at the specified index.
        
    Raises
    ------
    IndexError
        If the index is out of range for the entries in the toctree.
    """
    lines = read_rst_file(filepath)
    start, entry_start, entry_end = find_toctree_block(lines)
    entries = parse_toctree_entries(lines, entry_start, entry_end)
    
    if 0 < index <= len(entries):
        return os.path.basename(entries[index - 1])
    else:
        raise IndexError(f"Index {index} is out of range for toctree entries.")

# Examples:
# update_rst("examples.rst", action="add", new_entry="examples/new_example")
# update_rst("examples.rst", action="delete", target="examples/8fad")
# update_rst("examples.rst", action="insert", target="examples/env", new_entry="examples/new_between")

