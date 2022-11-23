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
    Generic method to save our files in similar manner
    :param text: The text that should be saved in the file
    :param filename: The name of the file we want to save
    :param path: The path to save the file to
    :param timestamp: timestamp of initialization of the pipeline
    :param add_timestamp: Boolean to add timestamp to the filename
    :param overwrite_check: Boolean to check if the file already exist.
    :return: None, only saved dedicated file.
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






