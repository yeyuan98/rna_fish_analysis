"""
    This script fetches the following information for each image of a given probe:
        mask path, fishdot path, sample name, pixel size x, pixel size y, pixel size z
"""


def summarize(sample_paths, output_file_path, wildcards, config):  # TODO
    """
        Generates summary entry for each image file of a given probe.
    :param sample_paths: list of paths to different sample dirs
    :param output_file_path: path to the specified output file
    :param wildcards: wildcards
    :param config: snakemake config
    """
    # 1st, get sample names of the given probe using snake_helper -> a list
    # 2nd, get image names for each given (sample, probe) combination and generate a pandas df
    # iterate on df rows, get each data
    # write to csv as defined by path
    pass


try:
    snakemake
except NameError:
    raise ReferenceError("Mainflow integration_summarize is only compatible with snakemake script directive.")
summarize(snakemake.input, snakemake.output[0], snakemake.wildcards, snakemake.config)
