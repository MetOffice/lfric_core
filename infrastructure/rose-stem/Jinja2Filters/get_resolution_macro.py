#!/usr/bin/env python
# -*- coding: utf-8 -*-
##############################################################################
# Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
# For further details please refer to the file LICENCE.original which you
# should have received as part of this distribution.
##############################################################################
'''
Implements a Jinja2 filter to strip the resolutions info.
'''
import ast
import re
import utilities


def get_resolution_macro(call):
    '''
    Takes a string return list of resolutions to use
        crun: Number of runs to do in the crun.
        ...
    @param [in] context Jinja2 instance to run macro against.
    @return List resulting from retrieving the resolutions.
    '''
    # pylint: disable=too-many-locals, too-many-branches

    if call.find('(') == -1:
        arguments = []
    else:
        macro_arguments = call[call.index('(')+1:call.rindex(')')]
        arguments = utilities.dictionary_from_arguments(macro_arguments)

    normal_arguments = [argument for argument in arguments
                        if argument.find('=') == -1]
    keyword_arguments = [argument for argument in arguments
                         if argument.find('=') != -1]

    argument_list = []
    for argument in normal_arguments:
        argument_list.append(ast.literal_eval(argument.strip()))

    argument_dictionary = {}
    for argument in keyword_arguments:
        key, value = re.split(' *= *', argument)
        argument_dictionary[key.strip()] = ast.literal_eval(value.strip())

    if len(argument_list) >= 2:

        app_name = argument_list[0]
        configuration = argument_list[1]
        app_key = app_name + '_' + configuration

        argument_dictionary = {}
        for argument in keyword_arguments:
            key, value = re.split(' *= *', argument)
            argument_dictionary[key.strip()] = ast.literal_eval(value.strip())

        # Return info about the resolutions arguments
        resource_list = []
        if 'resolutions' in argument_dictionary.keys():
            resource_list = argument_dictionary['resolutions']

        resource_support_meshes = ['']
        if 'support_meshes' in argument_dictionary.keys():
            resource_support_meshes = argument_dictionary['support_meshes']

        # We now consider the possibility that entries for the resolutions
        # can either just be the resolution (i.e. a string) or a
        # (resolution, timestep) pair.  For the latter we generate a
        # dictionary relating the resolution to the timestep.
        resource_dictionary = {}
        for entry in resource_list:
            if isinstance(entry, (list, tuple)):
                if entry[0] not in resource_dictionary.keys():
                    resource_dictionary[entry[0]] = set()
                if len(entry) > 1:
                    if isinstance(entry[1], (int, float)):
                        resource_dictionary[entry[0]].add(entry[1])
                    if isinstance(entry[1], (list, tuple)):
                        resource_dictionary[entry[0]].update(entry[1])

        return_value = (app_key,
                        resource_list,
                        resource_dictionary,
                        resource_support_meshes)

    else:

        return_value = None, None, None, None

    return return_value
