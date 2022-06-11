"""Helper functions useful for the snakemake path wrangling"""
from os.path import join
from os import listdir
from os import name as sysname


def samples_get(config, probe):
    """
        Query sample names given probe
    :param config: snakemake config
    :param probe: probe of interest
    :return: list of sample names, NOT path
    """
    samples = listdir(join(config["data_roots"]["input"], probe))
    samples = [f for f in samples if not f.startswith(".")]
    return samples


def output_workflow(workflow, config):
    """
        Generates workflow-specific paths, internally used by output_sample_dirs
    :param workflow: workflow specified
    :param config: snakemake config
    :return: a string of workflow-specific path
    """
    result = ""
    if workflow == "convert":
        params = "bfconvert"
        result = join("convert", params)
    elif workflow == "fishdot":
        params = "+".join([str(v) for v in config["fishdot"].values()])  # Fishdot params can be numerical.
        params = config["pipeline_light_train"]["fishdot"] + "++" + params
        result = join("fishdot", params)
    elif workflow == "segmentation":
        params = "+".join(config["segmentation"].values())
        params = config["pipeline_light_train"]["segmentation"] + "++" + params
        result = join("segmentation", params)
    else:
        raise ValueError("Unsupported workflow")
    return result


def output_workflow_all(config):
    """
        Based on snake_helper to generate a single string representing all workflow parameters
    :param config: snakemake config
    :return: string
    """
    result = output_workflow("convert", config)
    result += "&" + output_workflow("segmentation", config)
    result += "&" + output_workflow("fishdot", config)
    result = result.replace("/", "_")
    if sysname == 'nt':
        result = result.replace("\\", "_")
    return result


def output_sample_dirs(workflow, config, target):
    """
        Generates output directories for a given workflow; optionally for a given workflow and target
    :param workflow: workflow specified
    :param config: snakemake config
    :param target: optional, probe x sample joined path combination
    :return: a list of directory paths, either len = 1 (given target) or len = variable (no specified target).
    """
    result = []
    if target is None:
        #  Need to generate output directories for ALL probe x sample combination
        probes = config["probes"]
        sampless = [samples_get(config, probe) for probe in probes]
        for i in range(len(probes)):
            samples = sampless[i]
            probe = probes[i]
            for sample in samples:
                result.append(join("results", probe, sample, output_workflow(workflow, config)))
    else:
        #  Generate only for the given probe x sample combination
        result.append(join(config["data_roots"]["output"], target, output_workflow(config, workflow)))
    return result


#  (probe, sample) -> paths to RAW input images
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


def output_integration_all(stage, config):
    """
        Generate integration workflow file paths.
    :param stage: the integration stage, "summarize", "qc" or "plot"
    :param config: snakemake config
    """
    probes = config['probes']
    if stage == "summarize":
        samples = [join("results", probe, "integration", output_workflow_all(config), "samples.csv")
                   for probe in probes]
        dots = [join("results", probe, "integration", output_workflow_all(config), "dots.csv")
                   for probe in probes]
        results = samples + dots
    elif stage == "qc":
        #  TODO: QC Plots
        raise NotImplementedError("")
    elif stage == "plot":
        results = [join("results",probe,"integration",output_workflow_all(config),"countPlot.pdf")
                   for probe in probes]
    else:
        raise NotImplementedError("Unsupported integration stage")
    return results
