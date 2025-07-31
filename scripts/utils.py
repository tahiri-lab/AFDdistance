import os

import matplotlib.pyplot as plt


def ensure_directories_exist(file_path):
    # Get the directory path from the file path
    dir_path = os.path.dirname(file_path)

    # Check if the directory exists, if not, create it
    if not os.path.exists(dir_path):
        os.makedirs(dir_path)  # This creates all necessary intermediate directories

    return dir_path


def generate_unique_filename(file_path):
    # Separate the file name from the file extension
    base_name, extension = os.path.splitext(file_path)

    # Check if the file already exists
    if not os.path.exists(file_path):
        ensure_directories_exist(file_path)
        return file_path

    # If the file exists, start adding numbers to the filename
    counter = 1
    new_file_path = f"{base_name}_{counter}{extension}"

    # Keep increasing the counter until the new file name doesn't exist
    while os.path.exists(new_file_path):
        counter += 1
        new_file_path = f"{base_name}_{counter}{extension}"

    return new_file_path


def save_fig_safe(file_path):
    unique_file = generate_unique_filename(file_path)
    plt.savefig(unique_file)
