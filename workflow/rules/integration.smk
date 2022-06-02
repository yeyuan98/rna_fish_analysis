from snake_helper import *

# INPUT FUNCTION (probe) -> full paths to cell segmentation masks
def get_probe_samples(wildcards, config):
    """
        Get a list of major result directory paths (results/{probe}/{sample}) to different samples of a given probe
    :param wildcards: input function
    :param config: snakemake config
    :return: a list of paths to dirs
    """
    sample_names = samples_get(config, wildcards.probe)
    result_base = join("results", wildcards.probe)
    return [join(result_base, sn) for sn in sample_names]


def output_workflow_all(config):
    """
        Based on snake_helper to generate a single string representing all workflow parameters
    :param config: snakemake config
    :return: string
    """
    result =  "convert_" + output_workflow("convert", config)
    result += "_segmentation_" + output_workflow("segmentation", config)
    result += "_fishdot_" + output_workflow("fishdot", config)
    return result


rule integration_summarize:
    """integrates all images in a specific probe x sample to get the final fishdot list"""
    input:
        lambda wildcards: get_probe_samples(wildcards, config = config)
    output:
        join("results", "{probe}", "integration", output_workflow_all(config), "samples.csv")
    threads:
        config["resources"]["threads"]["integration"]
    resources:
        mem_mb=config["resources"]["mem_mb"]["integration"],
        time=config["resources"]["max_time"]["integration"]
    conda:
        "../envs/image.yaml"
    script:
        "../scripts/mainflow_integration_summarize.py"


rule integration_overlap:
    """performs segmentation vs fishdot overlapping to get fist dots in segmentation mask"""
    input:
        join("results","{probe}","integration",output_workflow_all(config),"samples.csv")
    output:
        join("results","{probe}","integration",output_workflow_all(config),"dots.csv")
    threads:
        config["resources"]["threads"]["integration"]
    resources:
        mem_mb=config["resources"]["mem_mb"]["integration"],
        time=config["resources"]["max_time"]["integration"]
    conda:
        "../envs/image.yaml"
    script:
        "../scripts/mainflow_integration_overlap.py"
