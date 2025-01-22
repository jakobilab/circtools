#!/usr/bin/env python3

import os
import sys
import subprocess

def process_paths(paths):
    modified_paths = []
    for path in paths:
        if os.path.isabs(path):
            # Prefix absolute path with "/host_os"
            modified_path = os.path.join("/host_os", path.lstrip("/"))
            modified_paths.append(modified_path)
        # not an (absolut) path, add unmodified
        else:
            modified_paths.append(path)
    return modified_paths

if __name__ == "__main__":
    # print(sys.argv)
    if len(sys.argv) < 2:
        sys.exit(1)

    user_paths = sys.argv[1:]
    modified_paths = process_paths(user_paths)

    # print("DEBUG: " + " ".join(modified_paths))
    #
    # subprocess.run(["pwd"])
    # subprocess.run(["ls", "-la"])
    # subprocess.run(["ls", "/host_os/home/tjakobi/repos/jakobilab/long_read_circRNA/test_fastq/"])

    modified_paths.insert(0,'circtools')

    subprocess.run(modified_paths)
