#!/usr/bin/env python
# -*- coding: utf-8 -*-
##############################################################################
# (c) Crown copyright 2018 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
'''
Library of helper functions for Jinja2 macros.
'''
_CLOSE_BRACKET = {')': '(', '}': '{', ']': '['}


def dictionary_from_arguments(string):
    '''
    Splits a string of comma separated values ignoring any commas appearing
    between brackets. i.e. ()[]{}

    string (str) To be split.
    Returns the split components as a list.
    '''
    items = []
    item = ''
    bracket_nest = []

    for character in string:
        if character == ',' and bracket_nest == []:
            items.append(item)
            item = ''
            continue

        item += character
        if character in _CLOSE_BRACKET.values():
            bracket_nest.append(character)
        elif character in _CLOSE_BRACKET.keys():
            if _CLOSE_BRACKET[character] != bracket_nest[-1]:
                raise Exception('Mismatched brackets in string: '
                                + string)
            _ = bracket_nest.pop(-1)

    if item:
        items.append(item)

    return items
