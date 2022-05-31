"""
    Script for the format convert workflow
"""

# General I/O
import os
import subprocess
from mainflow_helpers import *


def sample_convert(inputs, output_dir):
    """
        Given input raw image full paths, convert and write to output_dir.
    :param inputs: input raw image absolute paths, list
    :param output_dir: output directory absolute path, string
    :return: none
    """
    original_names = [os.path.split(inp)[1] for inp in inputs]
    output_names = [os.path.splitext(original)[0] for original in original_names]
    output_names = [outname + ".ome.tif" for outname in output_names]
    output_paths = [os.path.join(output_dir, outname) for outname in output_names]
    os.makedirs(output_dir, exist_ok=True)
    for i in range(len(inputs)):
        in_path = inputs[i]
        out_path = output_paths[i]
        print_current_time("Running conversion for " + in_path)
        run = subprocess.run(["bfconvert", "-no-upgrade", in_path, out_path], capture_output=True)
        print_check_run_info(run)
        print_current_time("Done conversion")


try:
    snakemake
except NameError:
    raise ReferenceError("Mainflow convert is only compatible with snakemake script directive.")
sample_convert(snakemake.input, snakemake.output[0])
