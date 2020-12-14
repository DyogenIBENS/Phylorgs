#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Basic implementation of some sed/awk commands.

The regex language is the one from Python re module.

Allows to provide compact commands on the fly, for example given to argparse
"""

import re

REGEX_INT = re.compile(r'\d+')


def parse_simple_subst(command: str):
    """Make a substitution function from a sed-like expression.
    The pattern can be any valid python re regex."""
    if command.startswith('s'):
        splitter = command[1]  # Usually '/'
    _, pattern, repl, modifiers = command.split(splitter)
    flags = False
    if 'Ii' in modifiers:
        flags |= re.I
    if 'mM' in modifiers:
        flags | re.M
    regex = re.compile(pattern, flags)
    match_counts = REGEX_INT.search(modifiers)
    count = 1
    if 'g' in modifiers:
        count = 0
    elif match_counts:
        count = int(match_counts.group(0))

    def substituter(string): return regex.sub(repl, string, count)
    return substituter

