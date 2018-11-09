#!/usr/bin/env python
# -*- coding: utf-8 -*-
##############################################################################
# Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
# For further details please refer to the file LICENCE.original which you
# should have received as part of this distribution.
##############################################################################
'''
Implements a Jinja2 filter to convert a dictionary into assignment strings.
'''


def dict_to_assign(context, indent):
    '''
    Takes a dictionary and returns a string of assigments "key = value".

    @param [in] context Dictionary
    @return String resulting from setting the environment.
    '''
    env_variables = []
    for key, value in context.items():
        if not isinstance(value, dict):
            env_variables.append('%s = %s' % (key, value))

    joining_str = '\n' + indent
    return joining_str.join(env_variables)
