"""
    Script for the fishdot workflow
"""

# General I/O
import os
import subprocess
import time
import image_base
from mainflow_helpers import *

# fishdot
import generators


def fishdot(inputs, output_dir, wildcards, config):
    """
        Given input single channel (fish signal) image full paths, run dot analysis and write to output_dir
        snakemake config params:
            data_roots/generated_scripts, tools/{airlocalize, airlocalize_template_ini, airlocalize_template_m},
            fishdot, threads_per_sample/fishdot, mem_mb_per_sample/fishdot
    :param inputs: FULL paths to single channel images (fish signal)
    :param output_dir: Output directory
    :param wildcards: wildcards involved; this is used to generate script folder path
    :param config: snakemake config
    :return: none
    """
    #  Note - AIRLOCALIZE by default output file names to be the same as inputs, but extension swapped to par4 and loc4.
    os.mkdir(output_dir)
    #  First, generate scripts
    params = "+".join([wildcards.probe, wildcards.sample, wildcards.params])
    params = wildcards.light_train + "++" + params
    script_dir = os.path.join(config["data_roots"]["generated_scripts"], "AIRLOCALIZE", params)
    output_mpath = os.path.join(script_dir, "run.m")
    os.makedirs(script_dir, exist_ok=True)  # For generated script, snakemake will NOT be involved, need to make dirs.
    inputs_abspath = [os.path.abspath(inp) for inp in inputs]  # get absolute path
    output_dir_abspath = os.path.abspath(output_dir)
    print_current_time(f"Generating MATLAB scripts for {params}")
    generators.AIRLOCALIZE_gen(config, inputs_abspath, output_mpath, output_dir_abspath)
    #  Next, call MATLAB for execution
    print_current_time(f"Running MATLAB for {params}")
    matlab_run_string = f"run('{output_mpath}'); exit;"
    run = subprocess.run(["matlab", "-nodisplay -nosplash -nodesktop", "-r", matlab_run_string])
    print_check_run_info(run)
    print_current_time("Done fishdot")


try:
    snakemake
except NameError:
    raise ReferenceError("Mainflow fishdot is only compatible with snakemake script directive.")
fishdot(snakemake.input, snakemake.output[0], snakemake.wildcards, snakemake.config)
