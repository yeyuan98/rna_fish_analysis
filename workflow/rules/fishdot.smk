from snake_helper import *


# INPUT FUNCTION (probe, sample, channel) -> full paths to single channel images
def get_singlechannel_sample(wildcards):
    """
        Get FULL paths to all single channel images for a given probe x sample x channel
    :param wildcards: input function
    :return: a list of FULL paths
    """
    probe = wildcards.probe
    sample = wildcards.sample
    train = wildcards.light_train
    singlechannel_outdir = checkpoints.sample_singlechannel.get(probe = probe, sample = sample, light_train = train).output[0]
    singlechannel_files = listdir(singlechannel_outdir)
    singlechannel_files = [f for f in singlechannel_files if not f.startswith(".")]
    return [join(singlechannel_outdir, sfile) for sfile in singlechannel_files]


rule fishdot:
    """singlechannel images -> dir. Python-end image processing is only used for fetching metadata of physical sizes."""
    input:
        get_singlechannel_sample
    output:
        directory(join("results", "{probe}", "{sample}", "fishdot", "{light_train}++{params}"))
    threads:
        config["resources"]["threads"]["fishdot"]
    resources:
        mem_mb=config["resources"]["mem_mb"]["fishdot"],
        time=config["resources"]["max_time"]["fishdot"]
    conda:
        "../envs/image.yaml"
    script:
        "../scripts/mainflow_fishdot.py"
