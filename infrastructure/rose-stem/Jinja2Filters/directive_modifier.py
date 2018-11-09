#!/usr/bin/env python2
# -*- coding: utf-8 -*-
##############################################################################
# Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
# For further details please refer to the file LICENCE.original which you
# should have received as part of this distribution.
##############################################################################
'''
Implements a Jinja 2 filter which inserts values into batch manager directive
strings.
'''

from __future__ import print_function

import math
import re


def directive_modifier(directive, cores, walltime, xios_nodes=0):
    '''
    Substitutes values into batch manager directives.
    '''
    def choose_replacement(name, arguments, xios_nodes):
        '''
        Converts a substitution name into a value.
        '''
        if name == 'nodes':
            node_size = int(arguments[0])
            middlebit = str(int(math.ceil(float(cores)
                                          / float(node_size))) + xios_nodes)
        elif name == 'cores':
            middlebit = str(cores)
        elif name == 'time_hhmmss':
            middlebit = str(walltime)
        else:
            raise Exception('Unrecognised function name "{}"'.format(name))
        return middlebit

    if isinstance(directive, str):
        directive = [directive]

    pattern = re.compile(r'<<(\S+?)\((\S*?)\)>>')

    plist = []
    for directive_i in directive:
        processed = directive_i
        for match in pattern.finditer(directive_i):
            name = match.group(1)
            arguments = match.group(2).split(',')
            match_string = '<<%s(%s)>>' % (match.group(1), match.group(2))
            processed = processed.replace(match_string,
                                          choose_replacement(name,
                                                             arguments,
                                                             xios_nodes))
        plist.append(processed)

    return plist


if __name__ == '__main__':
    DIRECTIVE = '-l=select=<<nodes(36)>>,walltime=<<time_hhmmss()>>'
    print(directive_modifier(DIRECTIVE, 122, '00:05:00'))
