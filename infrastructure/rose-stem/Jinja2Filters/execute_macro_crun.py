#!/usr/bin/env python
# -*- coding: utf-8 -*-
##############################################################################
# Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
# For further details please refer to the file LICENCE.original which you
# should have received as part of this distribution.
##############################################################################
'''
Implements a Jinja2 filter to run a macro specified by a string.
'''
import re
from jinja2 import contextfilter


@contextfilter
def execute_macro_crun(context, call):
    '''
    Takes a string and executes it as though it were a Jinja2 macro call, but
    overrides keyword do_crun to be True

    The call string has the syntax <macro name>([<argument>]...).

    Arguments can be either position or keyword based.

    @param [inout] context Jinja2 instance to run macro against.
    @param [in]    call    Invokation string.
    @return String resulting from calling the macro.
    '''
    if call.find('(') == -1:
        macro_name = call
        arguments = []
    else:
        macro_name = call[:call.index('(')]
        arguments = re.split(', *', call[call.index('(')+1:call.rindex(')')])

    normal_arguments = [argument for argument in arguments
                        if argument.find('=') == -1]
    keyword_arguments = [argument for argument in arguments
                         if argument.find('=') != -1]

    argument_list = []
    for argument in normal_arguments:
        if argument[0] == '"':
            argument_list.append(argument[1:-2])
        else:
            argument_list.append(argument)

    argument_dictionary = {}
    for argument in keyword_arguments:
        key, value = re.split(' *= *', argument)
        argument_dictionary[key] = value

    argument_dictionary['do_crun'] = True
    return context.vars[macro_name](*argument_list, **argument_dictionary)
