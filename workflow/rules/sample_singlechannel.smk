from os import listdir
from os.path import join


# INPUT FUNCTION (probe, sample) -> paths to CONVERTED images
def get_converted_sample(wildcards):
    """
        Get FULL paths to converted images for a given probe x sample
    :param wildcards: input function
    :return: a list of FULL paths
    """
    convert_outdir = checkpoints.sample_convert.get(probe = wildcards.probe, sample = wildcards.sample).output[0]
    convert_files = listdir(convert_outdir)
    convert_files = [f for f in convert_files if not f.startswith(".")]
    return [join(convert_outdir, cfile) for cfile in convert_files]

checkpoint sample_singlechannel:
    """converted images -> dir"""
    input:
        get_converted_sample
    output:
        directory(join("results", "{probe}", "{sample}", "single_channel", "{light_train}"))
    threads:
        config["resources"]["threads"]["sample_singlechannel"]
    resources:
        mem_mb=config["resources"]["mem_mb"]["sample_singlechannel"],
        time=config["resources"]["max_time"]["sample_singlechannel"]
    conda:
        "../envs/image.yaml"
    script:
        "../scripts/mainflow_singlechannel.py"
