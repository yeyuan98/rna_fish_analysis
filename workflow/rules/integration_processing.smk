from snake_helper import *

rule integration_countPlot:
    """
        getting the counts/cell ~ group plot
        this rule is responsible for generating all plots in the countPlots folder.
        merged.pdf -> count plot of all samples
        batch.qc.pdf -> count plot of all samples, faceted by batch
    """
    input:
        samples=join("results", "{probe}", "integration", output_workflow_all(config), "samples.csv"),
        plot=join("results", "{probe}", "integration", output_workflow_all(config), "plot.csv"),
        dots=join("results", "{probe}", "integration", output_workflow_all(config), "dots.csv")
    output:
        plot=join("results","{probe}","integration",output_workflow_all(config), "countPlots", "{plot_type}.pdf")
    threads:
        config["resources"]["threads"]["integration"]
    resources:
        mem_mb=config["resources"]["mem_mb"]["integration"],
        time=config["resources"]["max_time"]["integration"]
    conda:
        "../envs/plot.yaml"
    script:
        "../scripts/mainflow_integration_countPlot.R"
