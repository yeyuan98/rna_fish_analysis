"""
    Handy helper functions for workflow scripts
"""
import time


def print_current_time(message=None):
    """
        Prints current time. Used for logging.
    :param message: optional message to show
    :return: none
    """
    current_time = time.strftime('%Y-%m-%d %H:%M:%S %Z', time.localtime(time.time()))
    if message is None:
        print(current_time)
    else:
        print(current_time + " " + message)


def print_check_run_info(run):
    """
        Print and checks run result. Used for verifying subprocess status.
    :param run: subprocess.run object
    """
    print("STDOUT------")
    print(run.stdout)
    print("STDERR------")
    print(run.stderr)
    run.check_returncode()
