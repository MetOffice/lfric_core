#!/usr/bin/env python3
##############################################################################
# (c) Crown copyright 2022 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
"""
Unit test namelist reader generator.
"""
import io
from pathlib import Path
from subprocess import run
from textwrap import dedent

import pytest

import configurator.namelistdescription as description

PICKER_EXE = 'rose_picker'

HERE = Path(__file__).resolve().parent


class TestNamelistMeta():
    """
    Tests writing namelist reading source.
    """
    def test_module_write_empty(self):
        # pylint: disable=no-self-use
        """
        Writing an empty record.
        """
        output_file = io.StringIO()

        uut = description.NamelistDescription('test')
        with pytest.raises(description.NamelistDescriptionException):
            uut.write_module(output_file)

    def test_module_write_one_of_each(self, tmp_path: Path):
        # pylint: disable=no-self-use
        """
        Writing all possible field types.
        """
        uut = description.NamelistDescription('test')
        uut.add_value('vint', 'integer')
        uut.add_value('dint', 'integer', 'default')
        uut.add_value('sint', 'integer', 'short')
        uut.add_value('lint', 'integer', 'long')
        uut.add_value('dlog', 'logical', 'default')
        uut.add_value('vreal', 'real')
        uut.add_value('dreal', 'real', 'default')
        uut.add_value('sreal', 'real', 'single')
        uut.add_value('lreal', 'real', 'double')
        uut.add_value('treal', 'real', 'second')
        uut.add_string('vstr')
        uut.add_string('dstr', configure_string_length='default')
        uut.add_string('fstr', configure_string_length='filename')
        uut.add_enumeration('enum', enumerators=['one', 'two', 'three'])
        output_file = tmp_path / 'one_of_each_mod.f90'
        uut.write_module(output_file)

        expected_file = HERE / 'one_each_mod.f90'
        assert (expected_file.read_text(encoding='ascii')
                == output_file.read_text(encoding='ascii') + '\n')

    def test_module_write_growing(self, tmp_path: Path):
        # pylint: disable=no-self-use
        """
        Rewriting a growing record.
        """
        output_file = tmp_path / 'growing_mod.f90'

        uut = description.NamelistDescription('test')
        uut.add_value('foo', 'integer')
        uut.write_module(output_file)

        expected_file = HERE / 'first_growing_mod.f90'
        assert (output_file.read_text(encoding='ascii') + '\n'
                == expected_file.read_text(encoding='ascii'))

        uut.add_value('bar', 'real', 'default')
        uut.write_module(output_file)

        expected_file = HERE / 'second_growing_mod.f90'
        assert (output_file.read_text(encoding='ascii') + '\n'
                == expected_file.read_text(encoding='ascii'))

    def test_enumeration_only(self, tmp_path: Path):
        # pylint: disable=no-self-use
        """
        Writing only an enumeration field.
        """
        uut = description.NamelistDescription('enum')
        uut.add_enumeration('value', enumerators=['one', 'two', 'three'])
        output_file = tmp_path / 'enum_mod.f90'
        uut.write_module(output_file)

        expected_file = HERE / 'enum_only_mod.f90'
        assert (output_file.read_text(encoding='ascii') + '\n'
                == expected_file.read_text(encoding='ascii'))

    def test_more_than_one_enumeration(self, tmp_path: Path):
        # pylint: disable=no-self-use
        """
        Writing multiple enumeration fields.
        """
        uut = description.NamelistDescription('twoenum')
        uut.add_enumeration('first', enumerators=['one', 'two', 'three'])
        uut.add_enumeration('second', enumerators=['ay', 'bee', 'see'])
        output_file = tmp_path / 'multi_enum_mod.f90'
        uut.write_module(output_file)

        expected_file = HERE / 'two_enum_mod.f90'
        assert (output_file.read_text(encoding='ascii') + '\n'
                == expected_file.read_text(encoding='ascii'))

    def test_module_write_string(self, tmp_path: Path):
        # pylint: disable=no-self-use
        """
        Writing string value fields.
        """
        output_file = tmp_path / 'string_mod.f90'

        test_unit = description.NamelistDescription('mirth')
        test_unit.add_string('chuckle', 'default')
        test_unit.add_string('guffaw', 'default', '3')
        test_unit.add_string('hysterics', 'default', ':')
        test_unit.add_string('chortle', 'default', 'namelist:random=biggles')
        test_unit.write_module(output_file)

        expected_file = HERE / 'string_mod.f90'
        assert (output_file.read_text(encoding='ascii') + '\n'
                == expected_file.read_text(encoding='ascii'))

    def test_module_write_computed(self, tmp_path: Path):
        # pylint: disable=no-self-use
        """
        Writing calculated fields.
        """
        output_file = tmp_path / 'computed_mod.f90'

        uut = description.NamelistDescription('teapot')
        uut.add_value('foo', 'real', 'default')
        uut.add_value('fum', 'real', 'default')
        uut.add_computed('bar', 'real', 'foo ** 2', configure_kind='default')
        uut.write_module(output_file)

        expected_file = HERE / 'computed_mod.f90'
        assert (output_file.read_text(encoding='ascii') + '\n'
                == expected_file.read_text(encoding='ascii'))

    def test_module_write_constant(self, tmp_path: Path):
        # pylint: disable=no-self-use
        """
        Writing calculated fields containing constants.
        """
        output_file = tmp_path / 'constants_mod.f90'

        uut = description.NamelistDescription('cheese')
        uut.add_usage('FUDGE', module='constants_mod')
        uut.add_value('fred', 'real', 'default')
        uut.add_computed('wilma', 'real', 'fred * FUDGE',
                         configure_kind='default')
        uut.write_module(output_file)

        expected_file = HERE / 'constants_mod.f90'
        assert (output_file.read_text(encoding='ascii') + '\n'
                == expected_file.read_text(encoding='ascii'))

    def test_module_write_array(self, tmp_path: Path):
        # pylint: disable=no-self-use
        """
        Writing array fields.
        """
        output_file = tmp_path / 'array_module.f90'

        uut = description.NamelistDescription('aerial')
        uut.add_usage('esize', module='wibble_mod')
        uut.add_value('lsize', 'integer', 'default')
        uut.add_string('absolute', bounds='5')
        uut.add_value('inlist', 'integer', bounds='lsize')
        uut.add_value('outlist', 'real', bounds='esize')
        uut.add_value('unknown', 'integer', bounds=':')
        uut.write_module(output_file)

        expected_file = HERE / 'array_mod.f90'
        assert (output_file.read_text(encoding='ascii') + '\n'
                == expected_file.read_text(encoding='ascii'))


class TestNamelistConfigDescription():
    """
    Tests reading a configuration metadata file.
    """
    @staticmethod
    def _compile_dictionary(namelists):

        dictionary = {}

        for namelist in namelists:

            parameters = {}
            for parameter in namelist.get_parameters():

                parameters[parameter.name] = \
                    [parameter.fortran_type.intrinsic_type,
                     parameter.fortran_type.kind]
                if parameter.get_configure_type() == 'enumeration':
                    parameters[parameter.name].extend(
                        list(parameter.mapping.keys()))
                elif parameter.get_configure_type() == 'computed':
                    parameters[parameter.name].append(
                        parameter.computation)
                elif parameter.get_configure_type() == 'array':
                    if isinstance(parameter.bounds, list):
                        bounds = parameter.bounds[0]
                    else:
                        bounds = parameter.bounds
                    parameters[parameter.name].append('(' + bounds + ')')

            dictionary[namelist.get_namelist_name()] = parameters

        return dictionary

    def test_parser_good_file(self, tmp_path: Path):
        """
        Reading well formed file.
        """
        input_file = tmp_path / 'test.ini'
        input_file.write_text(dedent('''
            [namelist:fred]

            [namelist:fred=first_thing]
            type=character

            [namelist:fred=second]
            type=integer

            [namelist:fred=filename]
            type=character
            !string_length=filename

            [namelist:fred=choices]
            !enumeration=true
            values=foo, bar, baz, qux
            '''))

        picker_command = [PICKER_EXE, str(input_file)]
        _ = run(picker_command, cwd=input_file.parent, check=True)

        json_file = input_file.with_suffix('.json')

        uut = description.NamelistConfigDescription()
        result = self._compile_dictionary(uut.process_config(json_file))

        assert {'fred': {'choices':     ['integer', 'i_def',
                                         'foo', 'bar', 'baz', 'qux'],
                         'filename':    ['character', 'str_max_filename'],
                         'first_thing': ['character', 'str_def'],
                         'second':      ['integer', 'i_def']}} == result

    def test_only_enumeration(self, tmp_path: Path):
        """
        Reading enumerated field.
        """
        input_file = tmp_path / 'input.ini'
        input_file.write_text(dedent('''
            [namelist:barney]

            [namelist:barney=stuff]
            !enumeration=true
            values=one, two, three
            '''))

        picker_command = [PICKER_EXE, str(input_file)]
        _ = run(picker_command, cwd=input_file.parent, check=True)

        json_file = input_file.with_suffix('.json')

        uut = description.NamelistConfigDescription()
        result = self._compile_dictionary(uut.process_config(json_file))

        assert {'barney': {'stuff': ['integer', 'i_def',
                                     'one', 'two', 'three']}} == result

    def test_non_enumeration_no_type(self, tmp_path: Path):
        # pylint: disable=no-self-use
        """
        Reading choice field without type.
        """
        input_file = tmp_path / 'test.init'
        input_file.write_text(dedent('''
            [namelist:barney]

            [namelist:barney=stuff]
            values=one, two, three
            '''))

        picker_command = [PICKER_EXE, str(input_file)]
        _ = run(picker_command, cwd=input_file.parent, check=True)

        json_file = input_file.with_suffix('.json')

        uut = description.NamelistConfigDescription()
        with pytest.raises(description.NamelistDescriptionException):
            uut.process_config(json_file)

    def test_no_member_type(self, tmp_path: Path):
        # pylint: disable=no-self-use
        """
        Reading field without type.
        """
        input_file = tmp_path / 'input.init'
        input_file.write_text(dedent('''
            [namelist:santa]

            [namelist:santa=elf]
            length=:
            '''))

        picker_command = [PICKER_EXE, str(input_file)]
        _ = run(picker_command, cwd=input_file.parent, check=True)

        json_file = input_file.with_suffix('.json')

        uut = description.NamelistConfigDescription()
        with pytest.raises(description.NamelistDescriptionException):
            uut.process_config(json_file)

    def test_computed_fields(self):
        """
        Reading computed fields.
        """
        input_file = HERE / 'computed.ini'

        picker_command = [PICKER_EXE, str(input_file)]
        _ = run(picker_command, cwd=input_file.parent, check=True)

        json_file = input_file.with_suffix('.json')

        uut = description.NamelistConfigDescription()
        result = self._compile_dictionary(uut.process_config(json_file))

        assert {'teapot': {'foo': ['real', 'r_def'],
                           'fum': ['real', 'r_def'],
                           'bar': ['real', 'r_def', 'foo ** 2'],
                           'baz': ['real', 'r_def', 'PI * foo'],
                           'dosh': ['real', 'r_def',
                                    'milk + (foo ** 2) - (PI * fum)']}} \
            == result

    def test_constant_in_computed(self, tmp_path: Path):
        """
        Reading computed field with constants.
        """
        input_file = tmp_path / 'input.ini'
        input_file.write_text(dedent('''
            [namelist:cheese]

            [namelist:cheese=fred]
            type=real

            [!namelist:cheese=wilma]
            type=real
            expression=fred * source:constants_mod=FUDGE
            '''))

        picker_command = [PICKER_EXE, str(input_file)]
        _ = run(picker_command, cwd=input_file.parent, check=True)

        json_file = input_file.with_suffix('.json')

        uut = description.NamelistConfigDescription()
        result = self._compile_dictionary(uut.process_config(json_file))

        assert {'cheese': {'fred':  ['real', 'r_def'],
                'wilma': ['real', 'r_def', 'fred * FUDGE']}} == result

    def test_array_fields(self, tmp_path: Path):
        """
        Reading array fields.
        """
        input_file = tmp_path / 'input.ini'
        input_file.write_text(dedent('''
            [namelist:aerial]

            [namelist:aerial=fred]
            type=real

            [namelist:aerial=wilma]
            type=real
            length=:
            !bounds=source:constants_mod=FUDGE

            [namelist:aerial=betty]
            type=logical
            length=:
            !bounds=fred

            [namelist:aerial=dino]
            type=integer
            length=:
            !bounds=namelist:sugar=TABLET

            [namelist:aerial=bambam]
            type=integer
            length=:
            '''))

        picker_command = [PICKER_EXE, str(input_file)]
        _ = run(picker_command, cwd=input_file.parent, check=True)

        json_file = input_file.with_suffix('.json')

        uut = description.NamelistConfigDescription()
        result = self._compile_dictionary(
            uut.process_config(json_file))

        assert result == \
            {'aerial': {'bambam': ['integer', 'i_def', '(:)'],
                        'betty':  ['logical', 'l_def', '(fred)'],
                        'fred':   ['real', 'r_def'],
                        'wilma':  ['real', 'r_def', '(FUDGE)'],
                        'dino':   ['integer', 'i_def', '(TABLET)']}}
