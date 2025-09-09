# Author: ChatGPT 5

import re
import string
from typing import Dict, Optional

def _infer_pattern_from_format_spec(spec: str) -> str | None:
    """
    Infer regex from a Python format specifier string (simplified).
    Handles width and d/f type codes, e.g. 04d -> exactly 4 digits with leading zeros.
    """
    if not spec:
        return None

    # detect type (last char often type)
    type_char = spec[-1] if spec[-1].isalpha() else None

    # check for width
    m = re.match(r"([<>=^])?([0 ]?)(\d+)?", spec)
    align, fill, width = m.groups() if m else (None, None, None)

    if type_char in ("d", "b", "o", "x", "X"):
        if width:
            w = int(width)
            if spec.startswith("0"):  # zero-padded
                return rf"\d{{{w}}}"
            else:  # space-padded — just digits, but length = width
                return rf"\d{{1,{w}}}"
        return r"[+-]?\d+"

    if type_char in ("f", "e", "E", "g", "G"):
        # floats: allow width digits, optional decimal
        if width:
            return rf"[+-]?\d{{1,{width}}}(?:\.\d+)?"
        return r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?"

    return None

class FormatValidator:
    """
    Build a validator for a given .format-style template.
    Supports named (recommended) or positional fields.
    Repeated fields are enforced via backreferences.
    """

    def __init__(
        self,
        fmt: str,
        field_patterns: Optional[Dict[str, str]] = None,
        *,
        greedy: bool = False,
        flexible_ws: bool = False,
    ):
        """
        fmt: e.g. "Hello {name}, you are {age:d}"
        field_patterns: optional per-field regex (without anchors), e.g. {"name": r"[A-Za-z]+"}
        greedy: if True, fields use '.+' instead of '.+?' (rarely needed)
        flexible_ws: if True, runs of spaces in literals match any whitespace (\\s+)
        """
        self.fmt = fmt
        self.field_patterns = field_patterns or {}
        self.greedy = greedy
        self.flexible_ws = flexible_ws
        self._regex = self._compile()

    def _compile(self) -> re.Pattern:
        formatter = string.Formatter()
        parts = []
        seen_fields = set()
        default_pat = ".+" if self.greedy else ".+?"

        for literal, field_name, format_spec, conversion in formatter.parse(self.fmt):
            # Literal text
            if literal:
                lit = re.escape(literal)
                if self.flexible_ws:
                    # Turn escaped spaces into \s+ (keeps other spacing intact)
                    lit = re.sub(r"(?:\\ )+", r"\\s+", lit)
                parts.append(lit)

            if field_name is None:  # no field here
                continue

            # Normalize field name: allow positional {} or {0}
            if field_name == "":
                # auto-numbered positional fields aren't supported by str.format,
                # but handle defensively by assigning a synthetic name
                field_name = f"p{len(seen_fields)}"
            elif field_name.isdigit():
                field_name = f"p{field_name}"

            if field_name in seen_fields:
                # same field repeated -> enforce equality
                parts.append(f"(?P={field_name})")
                continue

            # Choose pattern: explicit > inferred-from-spec > default
            pat = (
                self.field_patterns.get(field_name)
                or _infer_pattern_from_format_spec(format_spec)
                or default_pat
            )
            parts.append(f"(?P<{field_name}>{pat})")
            seen_fields.add(field_name)

        regex = "^" + "".join(parts) + "$"
        return re.compile(regex)

    def fullmatch(self, s: str) -> bool:
        return self._regex.fullmatch(s) is not None

    def extract(self, s: str) -> Optional[Dict[str, str]]:
        m = self._regex.fullmatch(s)
        return m.groupdict() if m else None

    @property
    def pattern(self) -> str:
        return self._regex.pattern


# ---- Convenience functions ----

def validate_format(
    fmt: str,
    s: str,
    field_patterns: Optional[Dict[str, str]] = None,
    *,
    greedy: bool = False,
    flexible_ws: bool = False,
) -> bool:
    return FormatValidator(fmt, field_patterns, greedy=greedy, flexible_ws=flexible_ws).fullmatch(s)

def extract_from_format(
    fmt: str,
    s: str,
    field_patterns: Optional[Dict[str, str]] = None,
    *,
    greedy: bool = False,
    flexible_ws: bool = False,
) -> Optional[Dict[str, str]]:
    return FormatValidator(fmt, field_patterns, greedy=greedy, flexible_ws=flexible_ws).extract(s)
