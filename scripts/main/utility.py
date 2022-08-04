import random
import string
import os

letters = string.ascii_lowercase


class list_of_lists_container:
    def __init__(self, list_of_lists=[[""]]):
        self.list_of_lists = list_of_lists


    def set_list_of_lists(self, list_of_lists):
        self.list_of_lists = list_of_lists

    def get_list_of_lists(self):
        return self.list_of_lists

    def read_self(self, path_to_read):
        with open(path_to_read, 'r') as f:
            read_list_of_lists = f.read()
        self.list_of_lists = [item.split("#") for item in read_list_of_lists.split("$")]
        return self.list_of_lists

    def write_self(self, write_path):
        file_name = ''.join(random.choice(letters) for i in range(5))
        path_to_write = os.path.join(write_path, f'{file_name}.txt')
        to_write = "$".join(["#".join([(str(sub_item)) for sub_item in item]) for item in self.list_of_lists])
        with open(path_to_write,'w') as f:
            f.write(to_write)
        return path_to_write


class model_params_container:
    def __init__(self, model_params={}):
        self.model_params = model_params

    def set_model_params(self, model_params):
        self.model_params = model_params

    def get_model_params(self):
        return self.model_params

    def read_self(self, path_to_read):
        with open(path_to_read, 'r') as f:
            read_model_params = f.read()
        self.model_params = {item.split("#")[0]:item.split("#")[1].strip() for item in read_model_params.split("$")}
        return self.model_params

    def write_self(self, write_path):
        file_name = ''.join(random.choice(letters) for i in range(5))
        path_to_write = os.path.join(write_path, f'{file_name}.txt')
        with open(path_to_write,'w') as f:
            f.write("$".join([f'{key}#{val}' for key, val in self.model_params.items()]))
        return path_to_write
