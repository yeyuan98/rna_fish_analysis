from snake_helper import *
import json


#  INPUT FUNCTION (probe, sample) -> paths to RAW input images
def input_sample_images(wildcards, config):
    """
        Get FULL paths to all RAW input images given probe and sample
            Note - FULL path means relative to the workflow root folder
    :param wildcards: input function
    :param config: snakemake config
    :return: a list of paths
    """
    probe = wildcards.probe
    sample = wildcards.sample
    pathbase = join(config["data_roots"]["input"], probe, sample)
    images = listdir(pathbase)
    images = [f for f in images if not f.startswith(".")]
    return [join(pathbase, image) for image in images]


checkpoint sample_convert:
    """images -> dir"""
    input:
        lambda wildcards: input_sample_images(wildcards, config = config)
    output:
        directory(join("results", "{probe}","{sample}", output_workflow("convert", config)))
    threads:
        config["resources"]["threads"]["sample_convert"]
    resources:
        mem_mb=config["resources"]["mem_mb"]["sample_convert"],
        time=config["resources"]["max_time"]["sample_convert"]
    params:
        inputs_parsed=lambda wildcards, input: json.dumps(input),
        output_dir=lambda wildcards, output: output[0]
    container:
        "library://yeyuan/rna_fish/main:bftools_python"
    shell:
        """
        bfconvert_python '{params.inputs_parsed}' {params.output_dir}
        """
