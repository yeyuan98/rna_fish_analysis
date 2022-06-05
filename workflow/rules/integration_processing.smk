from snake_helper import *

rule integration_countPlot:
    """getting the counts/cell ~ group plot"""
    input:
        samples=join("results", "{probe}", "integration", output_workflow_all(config), "samples.csv"),
        plot=join("results", "{probe}", "integration", output_workflow_all(config), "plot.csv"),
        dots=join("results", "{probe}", "integration", output_workflow_all(config), "dots.csv")
    output:
        plot=join("results","{probe}","integration",output_workflow_all(config),"countPlot.pdf")
    threads:
        config["resources"]["threads"]["integration"]
    resources:
        mem_mb=config["resources"]["mem_mb"]["integration"],
        time=config["resources"]["max_time"]["integration"]
    conda:
        "../envs/plot.yaml"
    script:
        "../scripts/mainflow_integration_countPlot.R"