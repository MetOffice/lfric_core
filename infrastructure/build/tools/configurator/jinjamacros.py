#!/usr/bin/env python3
##############################################################################
# Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
# For further details please refer to the file LICENCE which you
# should have received as part of this distribution.
##############################################################################
"""
Jinja2 filters.
"""

from typing import List, Optional, Sequence


def decorate_macro(
    subject: Sequence[str],
    prefix: Optional[str] = None,
    postfix: Optional[str] = None,
) -> List[str]:
    """
    Bracket all input lines.

    :param subject: Lines of text to buttress.
    :param prefix: Put at the start of each line.
    :param postfix: Put at the end of each line.
    """
    result = list(subject)

    if prefix:
        result = [prefix + value for value in result]

    if postfix:
        result = [value + postfix for value in result]

    return result
