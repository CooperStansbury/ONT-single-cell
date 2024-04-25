import os
import sys
import subprocess
import shutil
import argparse




def check_directory(dir_path):
    """Checks if a directory exists.

    Args:
        dir_path (str): Path to the directory.

    Returns:
        str: The absolute path of the directory if it exists.

    Raises:
        NotADirectoryError: If the path exists but is not a directory.
        FileNotFoundError: If the directory does not exist.
    """
    if not os.path.exists(dir_path):
        raise FileNotFoundError(f"Directory does not exist: {dir_path}")

    if not os.path.isdir(dir_path):
        raise NotADirectoryError(f"Path is not a directory: {dir_path}")

    return os.path.abspath(dir_path)  # Return the absolute path


def ensure_dir_exists(dir_path, action='error'):
    """Checks if a directory exists, with options for behavior on existence.

    Args:
        dir_path (str): The path to the directory.
        action (str, optional): Action to take if the directory exists:
            * 'error': Raise a FileExistsError (default).
            * 'overwrite': Overwrite the existing directory.
            * 'ignore': Do nothing if the directory exists. 

    Raises:
        FileExistsError: If the directory exists and 'action' is set to 'error'.
    """

    if os.path.exists(dir_path):
        if action == 'error':
            raise FileExistsError(f"Directory already exists: {dir_path}")
        elif action == 'overwrite':
            shutil.rmtree(dir_path)
            os.makedirs(dir_path)
            print(f"Directory overwritten: {dir_path}")
        elif action == 'ignore':
            pass  # Do nothing
    else:
        os.makedirs(dir_path)
        print(f"Directory created: {dir_path}")



def demultiplex(fastq_dir, blaze_path):
    command = [
        blaze_path,
        
    ]
    print(command)
    # # Execute the command using subprocess.run
    # output = subprocess.run(" ".join(command), shell=True, capture_output=True)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Demultiplexing on the fly.")

    # Required arguments
    parser.add_argument("input_dir", help="Path to the input directory.")
    parser.add_argument("output_dir", help="Path to the output directory.")

    # Optional arguments
    parser.add_argument("--action", type=str, default='overwrite', 
                        help="Action for overwriting existing locations: 'ignore', 'error', 'overwrite'")
    
    args = parser.parse_args()

    # check input and output
    input_dir = check_directory(args.input_dir)
    output_dir = ensure_dir_exists(args.output_dir, 
                                   action=args.action)

    # # load the input directory
    # input_dir = sys.argv[1]
    # input_dir = check_directory(input_dir)

    # # check the subdir
    # fastq_dir = f"{input_dir}/fastq_pass/"
    # fastq_dir = check_directory(fastq_dir)

    # print(fastq_dir)
    