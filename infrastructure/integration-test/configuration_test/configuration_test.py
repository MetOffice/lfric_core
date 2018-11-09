#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-
##############################################################################
# (c) Crown copyright 2017 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
'''
Exercise namelist loading capabilities.
'''
from __future__ import print_function

import collections
import re
import sys

from testframework import MpiTest, TestEngine, TestFailed


##############################################################################
class configuration_test(MpiTest):
    '''
    Tests the ability to load namelists in an MPI environment.
    '''
    # pylint: disable=invalid-name

    _INJECT = 'configuration_test.nml'

    def __init__(self):
        super(configuration_test, self).__init__([sys.argv[1], self._INJECT],
                                                 processes=4)
        self._precision = 0.0005

    def test(self, return_code, out, err):
        # pylint: disable=too-many-locals

        if return_code != 0:
            print('Standard out: {out}'.format(out=out), file=sys.stderr)
            print('Standard error: {err}'.format(err=err), file=sys.stderr)
            message = 'Unexpected failure of test executable: {code}'
            raise TestFailed(message.format(code=return_code))

        expected = {'a_dim': 2,
                    'angle_deg': 7.998,
                    'angle_rad': 0.1396,
                    'an_enum': 'second',
                    'bounded_array_local_dim': [0.3, 0.4],
                    'bounded_array1_namelist_dim': [0.3, 0.4, 0.5],
                    'bounded_array2_namelist_dim': [0.5, 0.6, 0.7],
                    'bounded_array_source_dim': [0.1, 0.6, 0.7, 0.9],
                    'closed_array': [0.2, 0.3, 0.4],
                    'open_array': [1, 2, 3, 4, 5],
                    'some_string': 'chocolate teapot',
                    'whole_number': 13}

        variable_pattern = re.compile(r'(\d+)\s+(.+?)\s*:\s*(.+)\s*')
        list_pattern = re.compile(r'(\'.*?\'|".*?"|[^ ]+)')

        seen = collections.defaultdict(list)

        for process in range(0, 4):
            with open('result.{0}.txt'.format(process), 'r') as handle:
                contents = handle.readlines()

            for line in contents:
                match = variable_pattern.match(line)
                if match:
                    rank = int(match.group(1))
                    variable = match.group(2).strip()
                    seen[rank].append(variable)
                    value = match.group(3).strip()

                    if rank != process:
                        message = 'Found output for rank {0} in file from {1}'
                        raise TestFailed(message.format(rank, process))

                    list_match = list_pattern.findall(value)
                    if isinstance(list_match, list):
                        value = []
                    for item in list_match:
                        if item.startswith("'") or item.startswith('"'):
                            value.append(item[1:-1])
                        else:
                            value.append(item)

                    if variable not in expected:
                        message = 'Found unexpected variable "{variable}"' \
                                  ' in output'
                        raise TestFailed(message.format(variable=variable))

                    self._check(variable, expected[variable], value)

                additions = set(seen[0]) - set(expected.keys())
                if len(additions) > 0:
                    message = 'Addition variables found: {adds}'
                    raise TestFailed(message.format(adds=', '.join(additions)))

        if seen.keys() != [0, 1, 2, 3]:
            message = 'Incorrect ranks: ' + str(seen.keys())
            raise TestFailed(message)

        return 'One of each configuration type loaded'

    def _check(self, name, expected, found):
        # pylint: disable=consider-using-enumerate
        if not isinstance(expected, list):
            expected = [expected]
        if not isinstance(found, list):
            found = [found]

        if len(found) != len(expected):
            message = 'Expected variable "{variable}" to hold "{expect}"' \
                      ' but found "{found}", different number of values.'
            raise TestFailed(message.format(variable=name,
                                            expect=expected,
                                            found=found))

        for index in range(len(expected)):
            if isinstance(expected[index], float):
                maybe_zero = abs(float(found[index]) - expected[index])
                if maybe_zero > self._precision:
                    message = 'Expected index {index} of variable' \
                              ' "{variable}" to be within {precision} of' \
                              ' {expected} but found {found}'
                    raise TestFailed(message.format(variable=name,
                                                    index=index,
                                                    precision=self._precision,
                                                    expected=expected[index],
                                                    found=found[index]))
            else:  # expected[index] is not float
                if found[index] != str(expected[index]):
                    message = 'Expected index {index} of variable' \
                              '"{variable}" to hold "{expect}" but found' \
                              ' "{found}"'
                    raise TestFailed(message.format(variable=name,
                                                    index=index,
                                                    expect=expected[index],
                                                    found=found[index]))


##############################################################################
if __name__ == '__main__':
    TestEngine.run(configuration_test())
