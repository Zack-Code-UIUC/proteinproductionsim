import os
from datetime import datetime


def get_date_string():
    return datetime.now().strftime("%Y_%m_%d_%H_%M")


def make_directory_now(direc, name):
    """
    Create a directory based on the current working directory, the target directory and inputted name. A date and
    time string will be attached to the back of the name.

    Parameters:
    direc(string): target directory
    name(string) : directory name that need to be created.

    Return:
    string: the path of the created directory
    """
    # get current working directory
    cwd = os.getcwd()
    # construct the date and time string
    now = datetime.now()
    date_string = now.strftime("%Y_%m_%d_%H_%M")
    # make the directory name
    directory_name = name + '_ ' + date_string
    # construct the whole path
    path_directory = cwd + '/' + direc + '/' + directory_name
    try:
        os.makedirs(path_directory)
    except FileExistsError:
        pass
    return path_directory


def print_progress_bar(iteration, total):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
    """
    prefix = 'Progress:'
    suffix = 'Completed'
    fill = '█'
    length = 50
    decimals = 2

    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filled_length = int(length * iteration // total)
    bar = fill * filled_length + '-' * (length - filled_length)
    print(f'{prefix} |{bar}| {percent}% {suffix}', end='\r')

    if iteration == total:
        print()


def if_out_of_interval(num, low, high):
    if low < 0 and num >= high:
        return True
    elif high < 0 and num <= -low:
        return True
    if low <= num <= high:
        return False
    else:
        return True
