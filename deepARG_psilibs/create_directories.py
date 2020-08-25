#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-
import os
import sys


def create_dir(root_dir):
    subfolders = ["input", "model", "pssm_files",
                  "predictions", "database",
                  "prokka", "temp"]
    if not os.path.exists(root_dir):
        os.mkdir(root_dir)
    else:
        sys.stderr.write(
                'Cannot create folder "%s"! File/folder with '
                "the same name exists already.\n" % root_dir
            )
        sys.exit(2)
    for folder in subfolders:
        folder_dir = "{}/{}".format(root_dir, folder)
        if not os.path.exists(folder_dir):
            os.mkdir(folder_dir)
        else:
            sys.stderr.write(
                'Cannot create folder "%s"! File/folder with '
                "the same name exists already.\n" % folder_dir
            )
            sys.exit(2)
    sys.stdout.write(
        'Created folder "{}" and required subfolders.\n'.format(root_dir))

