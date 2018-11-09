#!/usr/bin/env python
# -*- coding: utf-8 -*-
##############################################################################
# Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
# For further details please refer to the file LICENCE.original which you
# should have received as part of this distribution.
##############################################################################
'''
Implements a Jinja2 filter which takes a number of strings and joins them with
the filter text.

e.g.

' foo ' | command_list('comm1', ['comm2, 'comm3'])
=> comm1 foo comm2 foo comm3
'''


def command_list(separator, *fragments):
    '''
    Takes an arbitrary collection of strings and joins them with the filter
    text into a single string.

    @param [in] separator Filter text string
    @param [in] fragments Multiple strings or lists of strings.
    '''
    commands = []
    for fragment in fragments:
        if hasattr(fragment, '__iter__'):
            commands.extend(fragment)
        else:
            commands.append(fragment)
    return separator.join(commands)
