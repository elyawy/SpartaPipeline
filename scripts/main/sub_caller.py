import subprocess
import tempfile
import os


class sub_caller:

    def __init__(self, process_path, output_path, length_dist):
        self.tempdir = tempfile.TemporaryDirectory()
        self.process_path = process_path
        self.output_path = output_path
        self.length_dist = length_dist


    def run_process(self, *args):
        tmpdirname = self.tempdir.name

        paths = []
        for arg in args:
            paths.append(arg.write_self(tmpdirname))

        p = subprocess.Popen(["python3", self.process_path, *paths, self.output_path, self.length_dist],
                            stdout=subprocess.PIPE, stdin=subprocess.PIPE
        )
        self.process = p
        print("started correction process, this may take a while")
    
    def get_result(self):
        result = self.process.communicate()[0].decode().strip()
        self.tempdir.cleanup()
        print("finished correction process")
        return result
