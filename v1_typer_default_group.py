# -*- coding: utf-8 -*-
"""
   adapted from click_default_group (see click_default_group description below)
   
   click_default_group
   ~~~~~~~~~~~~~~~~~~~
   Define a default subcommand by `default=True`:
   .. sourcecode:: python
      import click
      from click_default_group import DefaultGroup
      @click.group(cls=DefaultGroup, default_if_no_args=True)
      def cli():
          pass
      @cli.command(default=True)
      def foo():
          click.echo('foo')
      @cli.command()
      def bar():
          click.echo('bar')
   Then you can invoke that without explicit subcommand name:
   .. sourcecode:: console
      $ cli.py --help
      Usage: cli.py [OPTIONS] COMMAND [ARGS]...
      Options:
        --help    Show this message and exit.
      Command:
        foo*
        bar
      $ cli.py
      foo
      $ cli.py foo
      foo
      $ cli.py bar
      bar
"""
import warnings

import click
import typer

from typing import Any, Tuple


__all__ = ['DefaultGroup']
__version__ = '1.2.2'


class DefaultGroup(typer.Typer):
    """Invokes a subcommand marked with `default=True` if any subcommand not
    chosen.
    :param default_if_no_args: resolves to the default command if no arguments
                               passed.
    """

    def __init__(self, *args, **kwargs):
        # To resolve as the default command.
        if not kwargs.get('ignore_unknown_options', True):
            raise ValueError('Default group accepts unknown options')
        self.ignore_unknown_options = True
        self.default_cmd_name = kwargs.pop('default', None)
        self.default_if_no_args = kwargs.pop('default_if_no_args', False)
        self.commands = kwargs.pop('commands', {})
        self.params = kwargs.pop('params', {})
        print(self.commands)
        print(self.params)
        super(DefaultGroup, self).__init__(*args, **kwargs)
    
    def add_command(self, cmd, name = None):
        """Lifted from click src"""
        name = name or cmd.__name__
        if name is None:
            raise TypeError("Command has no name.")
        # _check_multicommand(self, name, cmd, register = True)
        self.commands[name] = cmd

    def set_default_command(self, command):
        """Sets a command function as the default command."""
        # cmd_name = command.name
        cmd_name = command.__name__
        self.add_command(command)
        # print(list(self.commands.keys()))
        self.default_cmd_name = cmd_name

    def parse_args(self, ctx, args):
        print("parse_args")
        if not args and self.default_if_no_args:
            args.insert(0, self.default_cmd_name)
        return super(DefaultGroup, self).parse_args(ctx, args)

    def get_command(self, ctx, cmd_name):
        print("we're in get_command")
        if cmd_name not in self.commands:
            print(f"get_command cmd_name: {cmd_name}")
            # No command name matched.
            ctx.arg0 = cmd_name
            cmd_name = self.default_cmd_name
        return super(DefaultGroup, self).get_command(ctx, cmd_name)

    def resolve_command(self, ctx, args):
        base = super(DefaultGroup, self)
        cmd_name, cmd, args = base.resolve_command(ctx, args)
        if hasattr(ctx, 'arg0'):
            args.insert(0, ctx.arg0)
            cmd_name = cmd.name
        return cmd_name, cmd, args

    def format_commands(self, ctx, formatter):
        formatter = DefaultCommandFormatter(self, formatter, mark='*')
        return super(DefaultGroup, self).format_commands(ctx, formatter)

    def command(self, *args, **kwargs):
        default = kwargs.pop('default', False)
        decorator = super(DefaultGroup, self).command(*args, **kwargs)
        if not default:
            return decorator
        warnings.warn('Use default param of DefaultGroup or '
                      'set_default_command() instead', DeprecationWarning)
        
        def _decorator(f):
            cmd = decorator(f)
            self.set_default_command(cmd)
            return cmd
        
        return _decorator
    
    def __call__(self, *args: Any, **kwargs: Any) -> Any:
        cmd = get_command(self)
        print(type(cmd))
        return cmd(*args, **kwargs)


class DefaultCommandFormatter(object):
    """Wraps a formatter to mark a default command."""

    def __init__(self, group, formatter, mark='*'):
        self.group = group
        self.formatter = formatter
        self.mark = mark

    def __getattr__(self, attr):
        return getattr(self.formatter, attr)

    def write_dl(self, rows, *args, **kwargs):
        rows_ = []
        for cmd_name, help in rows:
            if cmd_name == self.group.default_cmd_name:
                rows_.insert(0, (cmd_name + self.mark, help))
            else:
                rows_.append((cmd_name, help))
        return self.formatter.write_dl(rows_, *args, **kwargs)

def _check_multicommand(base_command, cmd_name, cmd, register = False):
    """Lifted from click src"""
    if not base_command.chain or not isinstance(cmd, MultiCommand):
        return
    if register:
        hint = (
            "It is not possible to add multi commands as children to"
            " another multi command that is in chain mode."
        )
    else:
        hint = (
            "Found a multi command as subcommand to a multi command"
            " that is in chain mode. This is not supported."
        )
    raise RuntimeError(
        f"{hint}. Command {base_command.__name__!r} is set to chain and"
        f" {cmd_name!r} was added as a subcommand but it in itself is a"
        f" multi command. ({cmd_name!r} is a {type(cmd).__name__}"
        f" within a chained {type(base_command).__name__} named"
        f" {base_command.__name__!r})."
    )

def get_command(typer_instance: DefaultGroup) -> click.Command:
    if typer_instance._add_completion:
        click_install_param, click_show_param = get_install_completion_arguments()
    if (
        typer_instance.registered_callback
        or typer_instance.info.callback
        or typer_instance.registered_groups
        or len(typer_instance.registered_commands) > 1
    ):
        print("ooh")
        # Create a Group
        click_command = typer.main.get_group(typer_instance)
        if typer_instance._add_completion:
            click_command.params.append(click_install_param)
            click_command.params.append(click_show_param)
        return click_command
    elif len(typer_instance.registered_commands) == 1:
        # Create a single Command
        click_command = get_command_from_info(typer_instance.registered_commands[0])
        if typer_instance._add_completion:
            click_command.params.append(click_install_param)
            click_command.params.append(click_show_param)
        return click_command
    assert False, "Could not get a command for this Typer instance"

def get_install_completion_arguments() -> Tuple[click.Parameter, click.Parameter]:
    install_param, show_param = typer.main.get_completion_inspect_parameters()
    click_install_param, _ = typer.main.get_click_param(install_param)
    click_show_param, _ = typer.main.get_click_param(show_param)
    return click_install_param, click_show_param
