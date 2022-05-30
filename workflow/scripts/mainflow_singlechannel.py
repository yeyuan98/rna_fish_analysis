"""
    Script for the single channel extraction workflow
"""
import os
import image_base


def sample_singlechannel(inputs, output_dir, light_train, config):
    """
        Given input converted image full paths, extract specified single channel and write to output_dir.
        snakemake config params:
            light_trains, threads_per_sample/sample_singlechannel, mem_mb_per_sample/sample_singlechannel
    :param inputs: list of paths to converted images
    :param output_dir: output directory
    :param light_train: which light train to extract
    :param config: snakemake config
    :return: none
    """
    original_names = [os.path.split(inp)[1] for inp in inputs]
    output_names = [original.replace(".ome.tif", ".tif") for original in original_names]
    output_paths = [os.path.join(output_dir, outname) for outname in output_names]
    os.mkdir(output_dir)
    try:
        train_query = config["light_trains"][light_train]
    except:
        raise KeyError("Specified light_train is undefined in config.")
    for i in range(len(output_paths)):
        in_path = inputs[i]
        out_path = output_paths[i]
        image_base.save_channel(in_path, out_path, train_query)


try:
    snakemake
except NameError:
    raise ReferenceError("Mainflow singlechannel is only compatible with snakemake script directive.")
sample_singlechannel(snakemake.input, snakemake.output[0], snakemake.wildcards.light_train, snakemake.config)
