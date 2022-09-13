########################################################################################
#  All generic functions that were used in different scripts
#  Made by Shane Pullens, Utrecht University - Theoretical bioinformatics.
#  Version 1.0
########################################################################################

from pathlib import Path
import os

def str2bool(v):
    return v.lower() in ("yes", "true", "t", "1", "True", "T", "Yes")


def write_file(text, filename, path, timestamp, add_timestamp=True, overwrite_check=False, verbose=False):
    """
    # TODO Add documentation
    :param text:
    :param filename:
    :param path:
    :param timestamp:
    :param add_timestamp:
    :param overwrite_check:
    :return:
    """
    if add_timestamp:  # Check if timestamp should be implemented
        filename = timestamp + '_' + filename  #+'.fasta' This was only here for testing right?

    # Check if the user already added an extension, skip adding another.
    filename = Path(filename)
    if filename.suffix not in ('.fasta', '.fa'):
        filename = str(filename) + '.fasta'

    absolute_path = Path(path / filename)

    # Check if the user enabled overwrite check. This checks if the user is going to overwrite a file and stops it.
    if overwrite_check:
        if os.path.isfile(absolute_path):
            print('>!! File: ' + str(filename) + ' already exists, aborted writing.')
            return None

    try:
        with open(absolute_path, 'w') as f:
            f.write(text)
        if verbose:
            print('>>> succesfully written '+filename)
    except:
        print('! An error occurred while trying to write '+filename)
    return None


