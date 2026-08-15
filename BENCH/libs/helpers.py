##########################################################
#  CROCO build system, under CeCILL-C
#  From Sébastien Valat (INRIA & LJK) - 2023
#  CROCO website : http://www.croco-ocean.org
##########################################################

##########################################################
"""Provide some helper function for example to run commands and patch files."""

##########################################################
import re
import os
import sys
import shutil
import subprocess
import tempfile
import traceback
from termcolor import cprint
from timeit import timeit
from .messaging import Messaging
from contextlib import contextmanager
from datetime import datetime


##########################################################
def display_run_error(command: str, outpout: str) -> None:
    dir = os.getcwd()
    keep_n_last_lines = 32
    truncated_output = "\n".join(outpout.split("\n")[-keep_n_last_lines:])
    cprint("-----------------------------------------------", "red")
    cprint(truncated_output, "dark_grey")
    cprint("-----------------------------------------------", "red")
    cprint(f"Error from command : {command}", "red")
    cprint(f"Error from workdir : {dir}", "red")
    cprint("-----------------------------------------------", "red")
    print("")


##########################################################
def print_exception(exception: Exception) -> None:
    dir = os.getcwd()
    cprint("-----------------------------------------------", "red")
    cprint("".join(traceback.format_exception(exception)), "dark_grey")
    cprint("-----------------------------------------------", "red")
    cprint(f"Error from command : {' '.join(sys.argv)}", "red")
    cprint(f"Error from workdir : {dir}", "red")
    cprint("-----------------------------------------------", "red")
    cprint(str(exception), "red")
    cprint("-----------------------------------------------", "red")
    print("")


##########################################################
def run_shell_command(command, logfilename=None, capture=True, show_on_error=True):
    """Print the command and run it. On failure it prints the output."""

    # print command
    Messaging.command(command)

    # capture or not capture
    if capture:
        # run
        result = subprocess.run(
            command,
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
        )
        if logfilename:
            with open(logfilename, "w") as log_file:
                log_file.write(result.stdout)
        # on failure
        if result.returncode != 0:
            if show_on_error:
                display_run_error(command, result.stdout)
            raise Exception("Fail to run command !")
    else:
        print("-----------------------------------------------")
        subprocess.run(command, shell=True, check=True)
        print("-----------------------------------------------")


def run_shell_command_time(command, logfilename=None, verbose: bool = False):
    """Print the command and run it by measuring time. On failure it prints the output."""
    # capture or not capture

    if verbose:
        return timeit(
            stmt=f"subprocess.run('{command}', check=True, shell=True)",
            setup="import subprocess",
            number=1,
        )
    else:
        if logfilename:
            log_fp = logfilename
        else:
            log_fp = tempfile.NamedTemporaryFile("w+").name

        try:
            return timeit(
                stmt=f"subprocess.run('({command}) 2>&1 > {log_fp}', check=True, shell=True)",
                setup="import subprocess",
                number=1,
            )
        except subprocess.CalledProcessError as e:
            # Handle subprocess-specific error
            display_run_error(command, log_fp)
            raise Exception(f"Command failed: {e}") from e
        except OSError as e:
            # Handle OS-related errors (e.g., command not found)
            display_run_error(command, log_fp)
            raise Exception(f"OS error occurred while running the command: {e}") from e
        except Exception as e:
            # Handle any other unforeseen exceptions
            display_run_error(command, log_fp)
            raise Exception(f"An error occurred while running the command: {e}") from e


def replace_in_file(in_path: str, out_path: str, pattern: str, replace_by: str):
    """Patch a file by replacing some patterns inside."""

    # print progress
    fname = os.path.basename(out_path)
    Messaging.step(f"Patching {fname}")

    # open & load
    with open(in_path, "r") as fp:
        content = fp.read()

    # patch
    content = content.replace(pattern, replace_by)

    # rewrite
    with open(out_path, "w+") as fp:
        fp.write(content)


def contains_one_of(seach_in: str, elements) -> bool:
    """Check is the string contains one of the given elements."""

    # if search a simple string
    if isinstance(elements, str):
        return elements in seach_in
    elif isinstance(elements, list):
        for element in elements:
            if element in seach_in:
                return True
        return False
    else:
        raise Exception("Unsupported type !")


def normalize_line(line):
    """Normalize a line by removing extra spaces and trimming."""
    if line is None:
        return ""
    return re.sub(r"\s+", " ", line).strip()


def apply_rule(fname, lines, rule):
    """
    Apply a single rule to a list of lines in memory.

    Parameters:
    - fname (str): The name of the file being patched.
    - lines (list): A list of strings (lines from the file).
    - rule (dict): A dictionary containing the rule to be applied. It can contain:
        - "what" (str, required): The line to match (after normalization) for certain operations.
        - "by" (str): Line used in 'replace' mode.
        - "insert" (str or list): Line(s) to insert in 'insert-before', 'insert-after', 'insert-after-before' modes.
        - "mode" (str, required): One of 'replace', 'insert-before', 'insert-after', 'insert-after-before', 'insert-at-begin', 'insert-at-end'.
        - "descr" (str, optional): Description of the rule.
        - "before" (str, required for 'insert-after-before'): The line before which to insert.
        - "after" (str, required for 'insert-after-before'): The line after which to insert.

    Returns:
    - list: The modified list of lines after applying the rule.
    """
    mode = rule.get("mode")
    what = normalize_line(rule.get("what", ""))
    descr = rule.get("descr", "No description provided")

    # Basic validations
    if not validate_rule(mode, what, rule, descr):
        return lines

    Messaging.step(f"Patching {fname} [ {descr} ]")

    # Determine the new line(s) to use based on the mode
    if mode == "replace":
        new_line = rule.get("by")
    else:
        new_line = rule.get("insert")

    # Handle rule application based on the mode
    if mode == "replace":
        return apply_replace_rule(lines, what, new_line, descr)
    elif mode == "insert-before":
        return apply_insert_before_rule(lines, what, new_line, descr)
    elif mode == "insert-after":
        return apply_insert_after_rule(lines, what, new_line, descr)
    elif mode == "insert-after-before":
        return apply_insert_after_before_rule(lines, rule, new_line, descr)
    elif mode == "insert-at-begin":
        return apply_insert_at_begin(lines, new_line, descr)
    elif mode == "insert-at-end":
        return apply_insert_at_end(lines, new_line, descr)

    print(f"Error: Unknown mode '{mode}' in rule '{descr}'. Skipping rule.")
    return lines


def validate_rule(mode, what, rule, descr):
    """Validate rule parameters."""
    if not mode:
        print(f"Error: No mode specified for rule '{descr}'. Skipping rule.")
        return False
    if mode not in [
        "replace",
        "insert-before",
        "insert-after",
        "insert-after-before",
        "insert-at-begin",
        "insert-at-end",
    ]:
        print(f"Error: Unknown mode '{mode}' in rule '{descr}'. Skipping rule.")
        return False
    if not what and mode not in ["insert-at-begin", "insert-at-end"]:
        print(
            f"Error: 'what' is required for mode '{mode}' in rule '{descr}'. Skipping rule."
        )
        return False
    if not rule.get("insert") and mode != "replace":
        print(
            f"Error: 'insert' is required for mode '{mode}' in rule '{descr}'. Skipping rule."
        )
        return False
    if mode == "insert-after-before":
        line_before = normalize_line(rule.get("before", ""))
        line_after = normalize_line(rule.get("after", ""))
        if not line_before or not line_after or not rule.get("insert"):
            print(
                f"Error: 'before', 'after', and 'insert' must be provided for 'insert-after-before' in rule '{descr}'. Skipping rule."
            )
            return False
    return True


def apply_replace_rule(lines, what, new_line, descr):
    """Apply 'replace' mode rule."""
    modified_lines = []
    inserted = False
    for line in lines:
        if normalize_line(line) == what:
            if isinstance(new_line, list):
                modified_lines.extend([nl + "\n" for nl in new_line])
            else:
                modified_lines.append(new_line + "\n")
            inserted = True
        else:
            modified_lines.append(line)
    if not inserted:
        print(f"Warning: No line matched '{what}' for replacement in rule '{descr}'.")
    return modified_lines


def apply_insert_before_rule(lines, what, new_line, descr):
    """Apply 'insert-before' mode rule."""
    modified_lines = []
    inserted = False
    for line in lines:
        if what in normalize_line(line) and not inserted:
            if isinstance(new_line, list):
                modified_lines.extend([nl + "\n" for nl in new_line])
            else:
                modified_lines.append(new_line + "\n")
            modified_lines.append(line)
            inserted = True
        else:
            modified_lines.append(line)
    if not inserted:
        print(f"Warning: No line matched '{what}' to insert before in rule '{descr}'.")
    return modified_lines


def apply_insert_after_rule(lines, what, new_line, descr):
    """Apply 'insert-after' mode rule."""
    modified_lines = []
    inserted = False
    for line in lines:
        if what in normalize_line(line) and not inserted:
            modified_lines.append(line)
            if isinstance(new_line, list):
                modified_lines.extend([nl + "\n" for nl in new_line])
            else:
                modified_lines.append(new_line + "\n")
            inserted = True
        else:
            modified_lines.append(line)
    if not inserted:
        print(f"Warning: No line matched '{what}' to insert after in rule '{descr}'.")
    return modified_lines


def apply_insert_after_before_rule(lines, rule, new_line, descr):
    """Apply 'insert-after-before' mode rule."""
    line_before = normalize_line(rule.get("before"))
    line_after = normalize_line(rule.get("after"))
    modified_lines = []
    inserted = False
    found_after = False
    for i, line in enumerate(lines):
        current_norm = normalize_line(line)
        modified_lines.append(line)

        if found_after and line_before in current_norm and not inserted:
            if isinstance(new_line, list):
                modified_lines[-1:-1] = [nl + "\n" for nl in new_line]
            else:
                modified_lines.insert(len(modified_lines) - 1, new_line + "\n")
            inserted = True

        if line_after in current_norm:
            found_after = True
    if not inserted:
        print(
            f"Warning: Could not find the sequence 'after' then 'before' for rule '{descr}'. No insertion performed."
        )
    return modified_lines


def apply_insert_at_begin(lines, new_line, descr):
    """Apply 'insert-at-begin' mode rule."""
    modified_lines = []
    if isinstance(new_line, list):
        modified_lines.extend([nl + "\n" for nl in new_line])
    else:
        modified_lines.append(new_line + "\n")
    modified_lines.extend(lines)
    return modified_lines


def apply_insert_at_end(lines, new_line, descr):
    """Apply 'insert-at-end' mode rule."""
    modified_lines = []
    modified_lines.extend(lines)
    if isinstance(new_line, list):
        modified_lines.extend([nl + "\n" for nl in new_line])
    else:
        modified_lines.append(new_line + "\n")
    return modified_lines


def patch_lines(file_path, rules):
    """
    Modify a file according to a list of rules.

    Parameters:
    - file_path: Path to the file to be modified.
    - rules: A list of dictionaries, each representing a rule.
    """
    # Backup the original file
    backup_path = f"{file_path}.backup"
    try:
        shutil.copyfile(file_path, backup_path)
        print(f"Backup created: {backup_path}")
    except FileNotFoundError:
        print(f"Error: File '{file_path}' not found. Aborting.")
        return
    except Exception as e:
        print(f"Error while creating backup: {e}")
        return

    # Read the original file
    try:
        with open(file_path, "r", encoding="utf-8") as f:
            lines = f.readlines()
    except Exception as e:
        print(f"Error reading the file: {e}")
        return

    # Apply each rule in-memory
    for rule in rules:
        lines = apply_rule(file_path, lines, rule)

    # Write the final result
    try:
        with open(file_path, "w", encoding="utf-8") as f:
            f.writelines(lines)
    #        print(f"Modifications completed successfully for {file_path}.")
    except Exception as e:
        print(f"Error writing the file: {e}")


@contextmanager
def move_in_dir(path):
    oldpwd = os.getcwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(oldpwd)


def apply_vars_in_str(value: str, vars: {}) -> str:
    # TODO: do not yet allow recursion but can ben done easly the way it is.

    # extract {case.value} patterns
    patterns = re.findall(r"\{[a-zA-Z0-9._-]+\}", value)

    # replace all
    if patterns:
        for pattern in patterns:
            # progress
            Messaging.step(f"Apply variable {pattern}...")

            # fetch value
            path = pattern.replace("{", "").replace("}", "").split(".")

            # start
            cursor = vars
            for element in path:
                cursor = cursor[element]

            # replace
            value = value.replace(pattern, cursor)

    # ok
    return value


def create_symlink(src_path, dest_path):
    """Helper function to create a symbolic link."""
    try:
        os.symlink(os.path.abspath(src_path), dest_path)
    except FileExistsError:
        print(f"The symbolic link {dest_path} already exists, link not created.")
    except OSError as e:
        print(f"Error creating symbolic link for {src_path}: {e}")


def process_symlink(src_path, dest_path):
    """Process and create a symlink for an existing symbolic link."""
    link_target = os.readlink(src_path)
    # Convert to absolute path if necessary
    link_target = os.path.join(os.path.dirname(src_path), link_target)
    link_target = os.path.abspath(link_target)

    if os.path.isfile(link_target):
        create_symlink(link_target, dest_path)


def copy_tree_with_absolute_symlinks(src, dest):
    """
    Copies a directory tree, recreating folders and replacing files with
    absolute symbolic links.
    Existing symbolic links in the source are recreated as symbolic links in
    the destination.


    :param src: Path to the source directory tree.
    :param dest: Path to the destination directory tree.
    """
    for root, dirs, files in os.walk(src, followlinks=True):
        # Compute the corresponding path in the destination
        relative_path = os.path.relpath(root, src)
        dest_dir = os.path.join(dest, relative_path)

        # Recreate folders in the destination
        os.makedirs(dest_dir, exist_ok=True)

        # Process directories and files in the current folder
        for name in dirs + files:
            src_path = os.path.join(root, name)
            dest_path = os.path.join(dest_dir, name)

            # Handle symbolic links and regular files
            if os.path.islink(src_path):
                process_symlink(src_path, dest_path)
            elif os.path.isfile(src_path):
                create_symlink(src_path, dest_path)
            # Directories are handled by os.makedirs above, so skip them


def extract_patch_by_partial_filename(patches, partial_filename):
    filtered_patches = {
        file: details for file, details in patches.items() if partial_filename in file
    }

    return filtered_patches


def suppress_patch(patches, patch_to_remove):
    target_file = patch_to_remove.get("file", "")
    target_after = patch_to_remove.get("after", "")
    if target_file in patches:
        patch_details = patches[target_file]
        # Vérifie si la ligne 'after' correspond
        if target_after in patch_details.get("after", ""):
            del patches[target_file]  # Supprime le patch du dictionnaire
            return True
    return False


def add_patch_with_list_support(patches, new_patch):
    file = new_patch["file"]
    patch_details = {
        "mode": new_patch.get("mode", "replace"),
        "after": new_patch.get("after", ""),
        "what": new_patch.get("what", ""),
        "by": new_patch.get("by", ""),
        "insert": new_patch.get("insert", ""),
        "descr": new_patch.get("descr", "No description provided."),
    }

    if file in patches:
        if isinstance(patches[file], list):
            patches[file].append(patch_details)
        else:
            patches[file] = [patches[file], patch_details]
    else:
        patches[file] = patch_details


def extract_elements_from_file(file_name, keyword, line_offset=1):
    """
    Extracts the elements of a line located a specified number of lines
    after a line containing a specific keyword in a file.

    :param file_name: The name of the file to read (str)
    :param keyword: The keyword to search for in a line (str)
    :param line_offset: The number of lines to skip after the keyword (int, default 1)
    :return: A list containing the elements of the target line
    """
    try:
        # Open the file and read all lines
        with open(file_name, "r") as file:
            lines = file.readlines()

        # Iterate through the lines to find the keyword
        for i, line in enumerate(lines):
            if line.strip().startswith(keyword):
                # Calculate the target line index
                target_line_index = i + line_offset
                if target_line_index < len(lines):  # Ensure the index is valid
                    target_line = lines[target_line_index].strip()
                    return target_line.split()

        # Return an empty list if the keyword is not found or no valid line exists
        return []
    except FileNotFoundError:
        print(f"Error: File '{file_name}' not found.")
        return []
    except Exception as e:
        print(f"An error occurred: {e}")
        return []


def delete_lines_from_file(file_name, keyword, line_offset=1, num_lines=1):
    """
    Deletes a specified number of lines from a file, starting from specific lines
    after a line identified by a keyword plus an optional offset (1-based offset).

    :param file_name: The name of the file to modify (str)
    :param keyword: The keyword to search for in the file (str)
    :param line_offset: The number of lines to skip after the keyword (1-based, int, default 1)
    :param num_lines: The number of lines to delete (int, default 1)
    """
    try:
        # Read all lines from the file
        with open(file_name, "r") as file:
            lines = file.readlines()

        # Find the index of the line containing the keyword
        for i, line in enumerate(lines):
            if keyword in line:
                # Calculate the first line to delete (adjust for 1-based offset)
                start_index = i + line_offset
                break
        else:
            print(f"Keyword '{keyword}' not found in the file.")
            return

        # Calculate the end index for deletion
        end_index = start_index + num_lines

        # Ensure the indices are within bounds
        if 0 <= start_index < len(lines) and 0 <= end_index <= len(lines):
            del lines[start_index:end_index]
        else:
            print("Error: Calculated indices are out of bounds.")
            return

        # Write the modified lines back to the file
        with open(file_name, "w") as file:
            file.writelines(lines)

        print(f"Deleted {num_lines} line(s) starting from line {start_index + 1}.")

    except FileNotFoundError:
        print(f"Error: File '{file_name}' not found.")
    except Exception as e:
        print(f"An error occurred: {e}")


def copy_and_replace(lst, index, new_value):
    """
    Creates a copy of the list and replaces the value at a specified index.

    :param lst: The original list (list)
    :param index: The index of the value to replace (int)
    :param new_value: The new value to insert (any type)
    :return: A new list with the replaced value
    """
    # Create a copy of the list
    new_list = lst.copy()

    # Replace the value at the given index
    if 0 <= index < len(new_list):  # Ensure index is within bounds
        new_list[index] = new_value
    else:
        print(f"Error: Index {index} is out of range.")

    return new_list


def delete_section_from_patch(patch, keyword):
    """
    Deletes a section from a patch dictionary starting from the keyword.

    :param patch: The dictionary containing patches (dict)
    :param keyword: The keyword to search for in the patch (str)
    """
    for file_path, operations in patch.items():
        # Filter out operations matching the keyword
        patch[file_path] = [
            op for op in operations if keyword not in op.get("what", "")
        ]

    print(f"Deleted sections containing keyword '{keyword}'.")


def parse_datetime(date_str):
    formats = ["%Y-%m-%d %H:%M:%S", "%Y/%m/%d %H:%M:%S"]

    for fmt in formats:
        try:
            return datetime.strptime(date_str, fmt)
        except ValueError:
            continue

    raise ValueError(f"Invalid date format: {date_str}")
