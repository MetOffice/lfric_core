#!/usr/bin/env python
# -*- coding: utf-8 -*-
##############################################################################
# (c) Crown copyright 2017 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
'''
Implements a Jinja2 filter which prints to standard error.
This is for debug purposes only.
'''

from __future__ import print_function

import sys


def debug_print(value):
    '''
    Sends the argument string to standard out for debug.
    '''
    print(value, file=sys.stderr)
