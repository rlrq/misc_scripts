# -*- coding: utf-8 -*-
"""
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
import inspect

import click
import typer

from typing import Any, Tuple

__all__ = ['DefaultGroup']
__version__ = '1.2.2'


class DefaultGroup(click.Group):
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
        super(DefaultGroup, self).__init__(*args, **kwargs)

    def set_default_command(self, command):
        """Sets a command function as the default command."""
        cmd_name = command.name
        self.add_command(command)
        self.default_cmd_name = cmd_name

    def parse_args(self, ctx, args):
        if not args and self.default_if_no_args:
            args.insert(0, self.default_cmd_name)
        return super(DefaultGroup, self).parse_args(ctx, args)

    def get_command(self, ctx, cmd_name):
        if cmd_name not in self.commands:
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

def get_command(typer_instance: DefaultGroup) -> click.Command:
    if typer_instance._add_completion:
        click_install_param, click_show_param = typer.main.get_install_completion_arguments()
    if (
        typer_instance.registered_callback
        or typer_instance.info.callback
        or typer_instance.registered_groups
        or len(typer_instance.registered_commands) > 1
    ):
        # Create a Group
        click_command = get_group(typer_instance)
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

def get_group(typer_instance: typer.Typer) -> click.Command:
    print("get grp")
    group = get_group_from_info(typer.models.TyperInfo(typer_instance))
    return group

def get_command_from_info(command_info: typer.models.CommandInfo) -> click.Command:
    assert command_info.callback, "A command must have a callback function"
    name = command_info.name or get_command_name(command_info.callback.__name__)
    use_help = command_info.help
    if use_help is None:
        use_help = inspect.getdoc(command_info.callback)
    else:
        use_help = inspect.cleandoc(use_help)
    (
        params,
        convertors,
        context_param_name,
    ) = typer.main.get_params_convertors_ctx_param_name_from_function(command_info.callback)
    cls = command_info.cls or TyperCommand
    command = cls(  # type: ignore
        name=name,
        context_settings=command_info.context_settings,
        callback=typer.main.get_callback(
            callback=command_info.callback,
            params=params,
            convertors=convertors,
            context_param_name=context_param_name,
        ),
        params=params,  # type: ignore
        help=use_help,
        epilog=command_info.epilog,
        short_help=command_info.short_help,
        options_metavar=command_info.options_metavar,
        add_help_option=command_info.add_help_option,
        no_args_is_help=command_info.no_args_is_help,
        hidden=command_info.hidden,
        deprecated=command_info.deprecated,
    )
    return command

def get_group_from_info(group_info: typer.models.TyperInfo) -> click.Command:
    assert (
        group_info.typer_instance
    ), "A Typer instance is needed to generate a Click Group"
    commands: Dict[str, click.Command] = {}
    for command_info in group_info.typer_instance.registered_commands:
        command = get_command_from_info(command_info=command_info)
        commands[command.name] = command
    for sub_group_info in group_info.typer_instance.registered_groups:
        sub_group = get_group_from_info(sub_group_info)
        commands[sub_group.name] = sub_group
    solved_info = typer.main.solve_typer_info_defaults(group_info)
    (
        params,
        convertors,
        context_param_name,
    ) = typer.main.get_params_convertors_ctx_param_name_from_function(solved_info.callback)
    # cls = solved_info.cls or click.Group
    cls = solved_info.cls or DefaultGroup
    group = cls(  # type: ignore
        name=solved_info.name or "",
        commands=commands,
        invoke_without_command=solved_info.invoke_without_command,
        no_args_is_help=solved_info.no_args_is_help,
        subcommand_metavar=solved_info.subcommand_metavar,
        chain=solved_info.chain,
        result_callback=solved_info.result_callback,
        context_settings=solved_info.context_settings,
        callback=typer.main.get_callback(
            callback=solved_info.callback,
            params=params,
            convertors=convertors,
            context_param_name=context_param_name,
        ),
        params=params,  # type: ignore
        help=solved_info.help,
        epilog=solved_info.epilog,
        short_help=solved_info.short_help,
        options_metavar=solved_info.options_metavar,
        add_help_option=solved_info.add_help_option,
        hidden=solved_info.hidden,
        deprecated=solved_info.deprecated,
    )
    print("grp commands:", group.commands)
    return group

