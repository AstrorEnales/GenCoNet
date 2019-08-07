import os
import shutil


def create_clean_directory(path: str):
    if path_is_existing_directory(path):
        delete_all_in_directory(path)
    else:
        os.mkdir(path)


def path_is_existing_directory(path: str) -> bool:
    return os.path.exists(path) and os.path.isdir(path)


def delete_all_in_directory(path: str):
    for child_name in os.listdir(path):
        child_path = os.path.join(path, child_name)
        try_delete_path(child_path)


def try_delete_path(path: str):
    try:
        delete_path(path)
    except Exception as e:
        print(e)


def delete_path(path: str):
    if os.path.isfile(path):
        os.unlink(path)
    elif os.path.isdir(path):
        shutil.rmtree(path)
