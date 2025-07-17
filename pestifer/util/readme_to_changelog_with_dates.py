# Author: ChatGPT 4o
"""
This script was written to do a one-time conversion of version notes in README.md to a CHANGELOG.md format.
It extracts version information from Git tags and PyPI releases, then formats it into a changelog.
"""

import subprocess
import requests
import re
from pathlib import Path

def get_git_tag_dates():
    """
    Get the creation dates of Git tags.
    This function retrieves the tags from the Git repository, sorts them by creation date,
    and returns a dictionary mapping tag names to their creation dates.
    The dates are formatted as 'YYYY-MM-DD'.
    """
    result = subprocess.run(
        ["git", "for-each-ref", "--sort=creatordate",
         "--format=%(refname:strip=2) %(creatordate:short)", "refs/tags"],
        capture_output=True, text=True
    )
    tag_date_map = {}
    for line in result.stdout.strip().splitlines():
        if " " in line:
            tag, date = line.strip().split(maxsplit=1)
            tag_date_map[tag.lstrip("v")] = date
    return tag_date_map

def get_pypi_release_dates(package_name: str) -> dict:
    """
    Get the release dates of a package from PyPI.
    This function queries the PyPI JSON API for the specified package,
    retrieves the release dates for each version, and returns a dictionary
    mapping version numbers to their release dates.
    The dates are formatted as 'YYYY-MM-DD'.
    """
    url = f"https://pypi.org/pypi/{package_name}/json"#!/usr/bin/env python3
    response = requests.get(url)
    response.raise_for_status()
    data = response.json()

    version_dates = {}
    for version, files in data.get("releases", {}).items():
        if files:
            upload_time = files[0].get("upload_time_iso_8601")
            if upload_time:
                version_dates[version] = upload_time[:10]  # YYYY-MM-DD
    return version_dates

def extract_changelog_from_readme(readme_path: Path, version_dates: dict) -> str:
    """
    Extract the changelog from the README.md file.
    This function reads the README.md file, identifies the release history section,
    and formats the version notes into a changelog format.
    It uses the provided version_dates dictionary to include release dates.
    The changelog is formatted with version headings and bullet points for each item.
    """
    with open(readme_path, "r") as f:
        lines = f.readlines()

    in_release_section = False
    changelog = ["# Changelog", "", "This project follows [Semantic Versioning](https://semver.org/) and documents changes below.", ""]
    current_version = None
    current_items = []

    for line in lines:
        if line.strip().lower().startswith("## release history"):
            in_release_section = True
            continue
        if in_release_section:
            version_match = re.match(r"^\* (\d+\.\d+(?:\.\d+)?[a-z0-9\-]*)", line.strip())
            if version_match:
                if current_version:
                    date_str = version_dates.get(current_version, "YYYY-MM-DD")
                    changelog.append(f"## [{current_version}] - {date_str}")
                    changelog.extend([f"- {item}" for item in current_items])
                    changelog.append("")
                current_version = version_match.group(1)
                current_items = []
            else:
                bullet_match = re.match(r"^\s*\*\s+(.*)", line)
                if bullet_match:
                    current_items.append(bullet_match.group(1).strip())
                elif line.strip() == "":
                    continue
                else:
                    if current_version:
                        date_str = version_dates.get(current_version, "YYYY-MM-DD")
                        changelog.append(f"## [{current_version}] - {date_str}")
                        changelog.extend([f"- {item}" for item in current_items])
                    break

    return "\n".join(changelog)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Convert README.md version notes to CHANGELOG.md using Git and PyPI")
    parser.add_argument("readme", type=Path, help="Path to README.md")
    parser.add_argument("--output", type=Path, default=Path("CHANGELOG.md"), help="Path to output CHANGELOG.md")
    args = parser.parse_args()

    tag_dates = get_git_tag_dates()
    pypi_dates = get_pypi_release_dates("pestifer")
    combined_dates = {**pypi_dates, **tag_dates}  # Git takes precedence

    changelog = extract_changelog_from_readme(args.readme, combined_dates)
    args.output.write_text(changelog)
    print(f"✅ CHANGELOG.md written to: {args.output.resolve()}")
