#!/usr/bin/env python3
##############################################################################
# (c) Crown copyright 2022 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
"""
Unit test feigner generator.
"""
from pathlib import Path

import configurator.namelistdescription as namelist
import configurator.namelistfeigner as feigner

HERE = Path(__file__).resolve().parent


###############################################################################
class TestFeigner:
    """
    Feigner generator unit tests.
    """
    def setUp(self):  # pylint: disable=invalid-name
        """
        Configures test framework.
        """
        self.maxDiff = None  # pylint: disable=invalid-name

    ###########################################################################
    def test_empty(self, tmp_path: Path):  # pylint: disable=no-self-use
        """
        Generation of feigner for empty namelist.
        """
        output_file = tmp_path / 'empty_mod.f90'
        uut = feigner.NamelistFeigner('empty_mod')
        uut.write_module(output_file)

        expected_file = HERE / 'empty_mod.f90'
        assert (expected_file.read_text(encoding='ascii')
                == output_file.read_text(encoding='ascii') + '\n')

    ###########################################################################
    def test_simple(self, tmp_path: Path):  # pylint: disable=no-self-use
        """
        Generation of feigner for simple value fields.
        """
        simple = namelist.NamelistDescription('simple')
        simple.add_value('foo', 'integer', 'default')
        simple.add_value('bar', 'real', 'double')
        simple.add_string('baz')
        simple.add_usage('qux', 'constants_mod')
        simple.add_value('fred', 'logical')

        output_file = tmp_path / 'feign.f90'
        uut = feigner.NamelistFeigner('simple_mod')
        uut.add_namelist([simple])
        uut.write_module(output_file)

        expected_file = HERE / 'simple_mod.f90'
        assert (expected_file.read_text(encoding='ascii')
                == output_file.read_text(encoding='ascii') + '\n')

    ###########################################################################
    def test_enumeration(self, tmp_path: Path):  # pylint: disable=no-self-use
        """
        Generation of feigner for enumerated fields.
        """
        enumable = namelist.NamelistDescription('enum')
        enumable.add_enumeration('thing', enumerators=['one', 'two'])

        output_file = tmp_path / 'enum_mod.f90'
        uut = feigner.NamelistFeigner('enumeration_mod')
        uut.add_namelist([enumable])
        uut.write_module(output_file)

        expected_file = HERE / 'enum_mod.f90'
        assert (expected_file.read_text(encoding='ascii')
                == output_file.read_text(encoding='ascii') + '\n')

    ###########################################################################
    def test_computed(self, tmp_path: Path):  # pylint: disable=no-self-use
        """
        Generation of feigner for computed fields.
        """
        simple = namelist.NamelistDescription('computed')
        simple.add_value('teapot', 'integer', 'default')
        simple.add_value('cheese', 'integer', 'default')
        simple.add_computed('biscuits', 'integer',
                            'teapot + cheese', 'default')

        output_file = tmp_path / 'computed_mod.f90'
        uut = feigner.NamelistFeigner('computed_mod')
        uut.add_namelist([simple])
        uut.write_module(output_file)

        expected_file = HERE / 'computed_mod.f90'
        assert (expected_file.read_text(encoding='ascii')
                == output_file.read_text(encoding='ascii') + '\n')

    ###########################################################################
    def test_everything(self, tmp_path: Path):  # pylint: disable=no-self-use
        """
        Generation of feigner for every kind of namelist value.
        """
        everything = namelist.NamelistDescription('everything')
        everything.add_string('cake', configure_string_length='filename')
        everything.add_enumeration('teapot', enumerators=['brown', 'steel'])
        everything.add_value('cheese', 'logical')
        everything.add_value('fish', 'real')
        everything.add_usage('wibble', 'constants_mod')
        everything.add_computed('yarn', 'real',
                                'fish * wibble / 180.0_r_def',
                                'default')
        everything.add_value('tail', 'integer', 'default')
        everything.add_value('school', 'integer', configure_kind='default',
                             bounds='2')
        everything.add_value('hanger', 'integer', configure_kind='default',
                             bounds='tail')
        everything.add_string('knife', configure_string_length='filename',
                              bounds='tail')

        output_file = tmp_path / 'everything_mod.f90'
        uut = feigner.NamelistFeigner('everything_mod')
        uut.add_namelist([everything])
        uut.write_module(output_file)

        expected_file = HERE / 'everything_mod.f90'
        assert (expected_file.read_text(encoding='ascii')
                == output_file.read_text(encoding='ascii') + '\n')

    ###########################################################################
    def test_multi_file(self, tmp_path: Path):  # pylint: disable=no-self-use
        """
        Generation of feigner for multiple namelists.
        """
        firstfile = namelist.NamelistDescription('first')
        firstfile.add_string('cake', configure_string_length='filename')
        firstfile.add_enumeration('teapot', enumerators=['brown', 'steel'])
        firstfile.add_value('cheese', 'logical')

        secondfile = namelist.NamelistDescription('second')
        secondfile.add_value('fish', 'real')
        secondfile.add_enumeration('yarn', enumerators=['fuzzy', 'colourful'])
        secondfile.add_value('tail', 'integer', 'default')

        output_file = tmp_path / 'multifile_mod.f90'
        uut = feigner.NamelistFeigner('multifile_mod')
        uut.add_namelist([firstfile, secondfile])
        uut.write_module(output_file)

        expected_file = HERE / 'multifile_mod.f90'
        assert (expected_file.read_text(encoding='ascii')
                == output_file.read_text(encoding='ascii') + '\n')
