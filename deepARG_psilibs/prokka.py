#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-

import subprocess
import sys
from Bio import SeqIO

def header_reducer(proj_root, fasta_file, out_prefix):
    outfile_path = "{}/input/{}.fa".format(proj_root,
                                           out_prefix)
    outfile = open(outfile_path, "w+")
    i = 0
    for record in SeqIO.parse(fasta_file, "fasta"):
        entries = record.description.split(" ")[1] + str(i)
        header = ">" + entries + "\n"
        outfile.write(header)
        outfile.write(str(record.seq) + "\n")
        i += 1
    outfile.close()
    return outfile_path



def prokka_translator(input_file, outdir, prefix):
    out_dir = "{}/prokka/{}".format(
                                outdir,
                                prefix)

    command = [
        "prokka",
        "--genus",
        "coxiella",
        "--usegenus",
        "--force",
        "--outdir",
        out_dir,
        "--prefix",
        prefix,
        input_file]
    subprocess.call(command)
    return "{}/{}.faa".format(out_dir, prefix)

