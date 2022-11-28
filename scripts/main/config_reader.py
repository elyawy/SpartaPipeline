

def parse_conf(path_to_conf_file):
    params = {}
    with open(path_to_conf_file,'r') as f:
        for line in f:
            name, value = line.strip().split(":")
            params[name.strip()] = value.strip()
    return params

