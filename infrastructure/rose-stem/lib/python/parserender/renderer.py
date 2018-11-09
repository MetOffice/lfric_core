#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-
##############################################################################
# Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
# For further details please refer to the file LICENCE.original which you
# should have received as part of this distribution.
##############################################################################
'''
Takes the results of parsing Cylc log files with parser.py and render them in
some attractive fashion.
'''
from __future__ import print_function

from abc import ABCMeta, abstractmethod
import datetime
from jinja2 import Environment, PackageLoader


##############################################################################
# Index rendering
##############################################################################
class IndexRenderer(object):
    '''
    Index page renderers inherit from this class.
    '''
    # pylint: disable=too-few-public-methods
    __MetaClass__ = ABCMeta

    @abstractmethod
    def render(self, title, file_list):
        '''
        Render to the specified file object.
        '''
        return ''


##############################################################################
class HtmlIndexRenderer(IndexRenderer):
    '''
    Renders a generic index page to an HTML document.
    '''
    # pylint: disable=too-few-public-methods
    def __init__(self):
        pass

    def render(self, title, file_list):
        variables = {'title': title,
                     'files': file_list}
        environment = Environment(loader=PackageLoader('parserender',
                                                       'templates'))
        template = environment.get_template('simpleindex.html')
        return template.render(variables)


##############################################################################
class SiteIndexRenderer(object):
    '''
    Index page renderers inherit from this class.
    '''
    # pylint: disable=too-few-public-methods
    __MetaClass__ = ABCMeta

    @abstractmethod
    def render(self, stream):
        '''
        Render to the specified file object.
        '''
        pass


##############################################################################
class HtmlSiteIndexRenderer(SiteIndexRenderer):
    '''
    Render the index page to an HTML document.
    '''
    # pylint: disable=too-few-public-methods

    def __init__(self, indexer, suite_url):
        self._indexer = indexer
        self._suite_url = suite_url

        if suite_url is not None:
            if self._suite_url.endswith('/'):
                self._suite_url = self._suite_url[:-1]

    def render(self, stream):
        threshold = datetime.datetime.utcnow() - datetime.timedelta(hours=48)
        variables = {'cronout': self._indexer.cronOut,
                     'crontimestamp': self._indexer.cronTimestamp,
                     'oldthreshold': threshold,
                     'suiteurl': self._suite_url,
                     'tree': self._indexer.getTree()}

        environment = Environment(extensions=['jinja2.ext.do',
                                              'jinja2.ext.loopcontrols'],
                                  loader=PackageLoader('parserender',
                                                       'templates'))
        template = environment.get_template('siteindex.html.jinja')
        print(template.render(variables), file=stream)


##############################################################################
# Compile output rendering
##############################################################################
class CompileRenderer(object):
    '''
    Compile event renderes inherit from this class.
    '''
    # pylint: disable=too-few-public-methods
    __MetaClass__ = ABCMeta

    @abstractmethod
    def render(self, ignore_codes, stream):
        '''
        Render to the specified file object.
        '''
        pass


##############################################################################
class HtmlCompileRenderer(CompileRenderer):
    '''
    Render the results of a compile to an HTML document.
    '''
    # pylint: disable=too-few-public-methods

    def __init__(self, context, status_parser, out_parser, error_parser):
        self._context = context
        self._status_parser = status_parser
        self._out_parser = out_parser
        self._error_parser = error_parser

    def render(self, ignore_codes, stream):
        if ignore_codes:
            event_buckets \
                = self._error_parser.get_events(ignore_codes=ignore_codes,
                                                group_by='filename')
        else:
            event_buckets = self._error_parser.get_events(group_by='filename')

        variables = {'context': self._context,
                     'compiler': self._out_parser.compiler,
                     'timestamp': self._status_parser.started,
                     'event_buckets': event_buckets}

        environment = Environment(loader=PackageLoader('parserender',
                                                       'templates'))
        template = environment.get_template('compilelog.html')
        print(template.render(variables), file=stream)


##############################################################################
# Run output rendering
##############################################################################
class RunRenderer(object):
    '''
    Run-time event renderers inherit from this class.
    '''
    # pylint: disable=too-few-public-methods
    __MetaClass__ = ABCMeta

    @abstractmethod
    def render(self, stream):
        '''
        Render to the specified file object.
        '''
        pass


##############################################################################
class HtmlRunRenderer(RunRenderer):
    '''
    Render the results of run-time checking to an HTML document.
    '''
    # pylint: disable=too-few-public-methods

    def __init__(self, context, out_parser, error_parser):
        self._context = context
        self._out_parser = out_parser
        self._error_parser = error_parser

    def render(self, stream):
        variables = {'context': self._context,
                     'compiler': self._out_parser.compiler,
                     'timestamp': self._out_parser.started,
                     'events': self._error_parser.get_events()}

        environment = Environment(loader=PackageLoader('parserender',
                                                       'templates'))
        template = environment.get_template('runlog.html')
        print(template.render(variables), file=stream)
