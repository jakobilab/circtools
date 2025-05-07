#!/usr/bin/env python3

import os
import sys
import subprocess
import shlex

def process_paths(paths):
    modified_paths_int = []
    for path in paths:
        if os.path.isabs(path):
            # Prefix absolute path with "/host_os"
            modified_path = os.path.join("/host_os", path.lstrip("/"))
            modified_paths_int.append(modified_path)
        # not an (absolut) path, add unmodified
        else:
            modified_paths_int.append(path)
    return modified_paths_int

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("✅ Patched version: no arguments provided.")
        sys.exit(1)

    user_paths = sys.argv[1:]
    modified_paths = process_paths(user_paths)

    args_str = ' '.join([shlex.quote(arg) for arg in modified_paths])
    command = f'. /circtools/bin/activate && circtools {args_str}'

    print(f"✅ Running: {command}")
    subprocess.run(command, shell=True)
