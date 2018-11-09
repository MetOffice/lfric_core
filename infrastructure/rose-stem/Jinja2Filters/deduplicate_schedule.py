#!/usr/bin/env python
# -*- coding: utf-8 -*-
##############################################################################
# Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
# For further details please refer to the file LICENCE.original which you
# should have received as part of this distribution.
##############################################################################
'''
Implements a Jinja2 filter which removes duplicates from a Cylc task schedule.
'''
import re

_LINE_TEMPLATE = '{indent}{prerequisite} => {result}'
_GROUP_HEAD_TEMPLATE = ['{indent1}{cyclegroup}',
                        '{indent2}graph="""']
_GROUP_FOOT_TEMPLATE = ['{indent2}"""']

_CYLC_PATTERN = re.compile(r'\s*(\[\[\[\S*\]\]\])(.*?)(?=\[\[\[)', re.S)
_DEPENDENCY_PATTERN = re.compile(r'\s*(\S+)\s*=>\s*(\S+)')


def deduplicate_schedule(schedule):
    '''
    Takes a Cylc task schedule and removes duplicates.  It also combines
    scheduling dependencies into groups with the same cycling pattern.

    For example, if schedule is:
       a => b
       e => f
       [[[R1]]]
           graph = """
           a => b
           c => d
           """
       [[[R3/P1/3]]]
           graph = """
           a[-P1] => a
           """
     then the returned schedule would be:
        [[[R1]]]
           graph = """
              a => b
              c => d
              e => f
              """
        [[[R3/P1/3]]]
           graph = """
              a[-P1] => a
              """
    i.e. tasks at the beginning of the schedule that aren't explicitly
         associated with a cycle period are put into the [[[R1]]] group of
         tasks.

    @param [in] schedule String containing Cylc task schedule.
    @return Deduplicated version of schedule.
    '''
    schedule_dict = {}
    for match in _CYLC_PATTERN.finditer('[[[]]]' + schedule + '[[['):
        key, value = match.group(1), match.group(2)
        if key == '[[[]]]':
            key = '[[[R1]]]'
        if key in schedule_dict.keys():
            schedule_dict[key] = schedule_dict[key] + value
        else:
            schedule_dict[key] = value

    new_schedule = []
    for key, value in schedule_dict.items():
        dependency_set = set()
        for match in _DEPENDENCY_PATTERN.finditer(value):
            dependency_set.add((match.group(1), match.group(2)))

        schedule_dict[key] = [_LINE_TEMPLATE.format(indent=' '*12,
                                                    prerequisite=first,
                                                    result=second)
                              for first, second in dependency_set]
        new_schedule.extend([line.format(indent1=' '*8,
                                         indent2=' '*12,
                                         cyclegroup=key)
                             for line in _GROUP_HEAD_TEMPLATE])
        new_schedule.extend(schedule_dict[key])
        new_schedule.extend([line.format(indent1=' '*8,
                                         indent2=' '*12,
                                         cyclegroup=key)
                             for line in _GROUP_FOOT_TEMPLATE])

    return '\n'.join(new_schedule)
