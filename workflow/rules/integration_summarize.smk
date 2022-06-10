from snake_helper import *

"""
    Generates summary tables for later integration steps.
    Rule integration_summarize is required for all other integration rules.
"""

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


rule integration_summarize:
    """integrates all images in a specific probe x sample to get the final fishdot list"""
    input:
        lambda wildcards: get_probe_samples(wildcards, config = config)
    output:
        samples=join("results", "{probe}", "integration", output_workflow_all(config), "samples.csv"),
        plot=join("results", "{probe}", "integration", output_workflow_all(config), "plot.csv"),
        plot_config=join("results", "{probe}", "integration", output_workflow_all(config), "plot.yaml"),
        visualization=directory(join("results", "{probe}", "integration", output_workflow_all(config), "visualization"))
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
