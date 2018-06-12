#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-
##############################################################################
# Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
# For further details please refer to the file LICENCE.original which you
# should have received as part of this distribution.
##############################################################################
'''
The Fortran logging module terminates on error. This cannot be tested by the
unit testing framework as it terminates the unit tests as well.
'''

from __future__ import print_function

import datetime
import re

from testframework import LFRicLoggingTest, MpiTest, TestEngine, TestFailed

##############################################################################
class log_mod_error_serial_test( MpiTest ):
  '''
  Tests that logging an error terminates execution when run serially.
  '''
  def __init__( self ):
    super(log_mod_error_serial_test, self).__init__( processes=1 )

    self._minimumTimestamp = datetime.datetime.utcnow()

  def test( self, returncode, out, err ):
    expectedLevel = 'ERROR'
    expectedMessage = ' An error was logged.'

    if returncode == 0:
      raise TestFailed( 'Logging an error did not cause termination to end' )

    if out != '':
      message = 'Expected no output on standard out but found: {out}' \
                + '\nStandard error: {err}'
      raise TestFailed( message.format( out=out, err=err ) )

    try:
      timestampString, level, report = err.split( ':', 3 )
      timestampWithoutTimezone = timestampString[:-5]

      timestamp = datetime.datetime.strptime( timestampWithoutTimezone, \
                                              '%Y%m%d%H%M%S.%f' )
    except Exception as ex:
      raise TestFailed( 'Unexpected log message: {}'.format( err ) )

    if timestamp < self._minimumTimestamp:
      message = 'Expected a timestamp after {} but read {}'
      raise TestFailed( message.format( minimumTimestamp, timestamp ) )

    if level != expectedLevel:
      message = 'Expected "{}" but read "{}"'
      raise TestFailed( message.format( expectedLevel, level ) )

    # We only check the first line as compilers tend to print the return code
    # as well. This will remain true until we can use Fortran 2008 and
    # "stop error".
    #
    first, newline, rest = report.partition( '\n' )
    if first != expectedMessage:
      message = 'Expected "{}" but read "{}"'
      raise TestFailed( message.format( expectedMessage, first ) )

    message = 'Logging an error caused exit as expected with code {code}'
    return message.format( code=returncode )

##############################################################################
class log_mod_error_parallel_test( LFRicLoggingTest ):
  '''
  Tests that logging an error terminates execution when run in parallel.
  '''
  def __init__( self ):
    super(log_mod_error_parallel_test, self).__init__( processes=2 )

    self._minimumTimestamp = datetime.datetime.utcnow()
    self._linePattern = re.compile( r'(\d\d\d\d)(\d\d)(\d\d) (\d\d)(\d\d)(\d\d)\.(\d\d\d) (\S+)\s+PET(\d+)\s+(.*)' )

  def test( self, returncode, out, err ):
    expectedLevel = 'ERROR'
    expectedMessage = 'An error was logged.'

    if returncode == 0:
      raise TestFailed( 'Logging an error did not cause termination to end' )

    if out != '' or err != '':
      message = 'Expected no output on standard out or standard error:\n' \
                 + 'Standard out: {out}\nStandard error: {err}'
      raise TestFailed( message.format( out=out, err=err ) )

    # We remove the first line as it will be a spin-up message.
    petLog = self.getLFRicLoggingLog()
    petLog = '\n'.join( petLog.splitlines()[1:] )

    match = self._linePattern.match( petLog )
    if match:
      try:
        timestamp = datetime.datetime( int(match.group(1)), # Year
                                       int(match.group(2)), # Month
                                       int(match.group(3)), # Day
                                       int(match.group(4)), # Hour
                                       int(match.group(5)), # Minute
                                       int(match.group(6)), # Second
                                    int(match.group(7)) * 1000, # Microseconds
                                       None ) # Timezone
      except Exception as ex:
        raise TestFailed( 'Bad timestamp format: {}'.format( petLog ) )
      level = match.group(8)
      report = match.group(10)
    else:
      raise TestFailed( 'Unexpected log message: {}'.format( petLog ) )

    if timestamp < self._minimumTimestamp:
      message = 'Expected a timestamp after {} but read {}'
      raise TestFailed( message.format( minimumTimestamp, timestamp ) )

    if level != expectedLevel:
      message = 'Expected "{}" but read "{}"'
      raise TestFailed( message.format( expectedLevel, level ) )

    # We only check the first line as compilers tend to print the return code
    # as well. This will remain true until we can use Fortran 2008 and
    # "stop error".
    #
    first, newline, rest = report.partition( '\n' )
    if first != expectedMessage:
      message = 'Expected "{}" but read "{}"'
      raise TestFailed( message.format( expectedMessage, first ) )

    message = 'Logging an error caused exit as expected with code {code}'
    return message.format( code=returncode )

##############################################################################
if __name__ == '__main__':
  TestEngine.run( log_mod_error_serial_test() )
  TestEngine.run( log_mod_error_parallel_test() )
