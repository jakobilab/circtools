#!/usr/bin/env python3

import os
import sys
import subprocess
import shlex

def process_paths(paths):
    modified_paths_int = []
    for path in paths:
        if os.path.isabs(path) or (path.startswith("@") and os.path.isabs(path.replace("@",""))):
            # Prefix absolute path with "/host_os"
            modified_path = os.path.join("/host_os", path.replace("@","").lstrip("/"))

            # check if macOS system
            # is yes, we have to fix the absolute path and insert /host_mnt/
            # i.e. /host_os/Users/tjakobi/tmp becomes
            #      /host_os/host_mnt/Users/tjakobi/tmp
            if os.path.isdir("/host_os/host_mnt/Users/"):

                # this is only necessary for Users and Volume paths
                modified_path = modified_path.replace("/host_os/Users/", "/host_os/host_mnt/Users/")
                modified_path = modified_path.replace("/host_os/Volumes/", "/host_os/host_mnt/Volumes/")

            if path.startswith("@"):
                modified_paths_int.append("@"+modified_path)
            else:
                modified_paths_int.append(modified_path)

        # not an (absolut) path, add unmodified
        else:
            modified_paths_int.append(path)

    return modified_paths_int

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("No arguments provided.")
        sys.exit(1)

    user_paths = sys.argv[1:]
    modified_paths = process_paths(user_paths)

    args_str = ' '.join([shlex.quote(arg) for arg in modified_paths])
    command = f'source /circtools/bin/activate && cd /host_rel/ && circtools {args_str}'

    subprocess.run(command, shell=True, executable="/bin/bash")
