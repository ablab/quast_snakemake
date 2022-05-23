import sys
from datetime import datetime


def print_error(msg, indent='', exit_code=None):
    sys.stderr.write(indent + 'ERROR! ' + msg + '\n')
    if exit_code:
        exit(exit_code)


def print_warning(msg, indent=''):
    sys.stderr.write(indent + 'WARNING! ' + msg + '\n')


def print_notice(msg, indent=''):
    sys.stderr.write(indent + 'NOTICE: ' + msg + '\n')


def print_info(msg, indent=''):
    sys.stdout.write(indent + 'INFO: ' + msg + '\n')


def print_debug(msg, indent=''):
    sys.stderr.write(indent + 'DEBUG: ' + msg + '\n')


def print_command_line(args, indent = '', wrap_after = 80, only_if_debug=True):
    text = ''
    line = indent

    for i, arg in enumerate(args):
        if ' ' in arg or '\t' in arg:
            args[i] = "'" + arg + "'"

        line += arg

        if i == len(args) - 1:
            text += line

        elif wrap_after is not None and len(line) > wrap_after:
            text += line + ' \\\n'
            line = ' ' * len(indent)

        else: line += ' '

    sys.stdout.write(text)


def print_timestamp(msg=''):
    now = datetime.now()
    current_time = now.strftime("%Y-%m-%d %H:%M:%S")
    sys.stdout.write('')
    sys.stdout.write(msg + current_time + '\n')

