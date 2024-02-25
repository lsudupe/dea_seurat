
rule process_data:
    input:
        lambda wildcards: f"/Users/medinils/Desktop/dea_seurat/test_data/{data_paths[wildcards.analysis]}"
    output:
        "results/{analysis}/processed_data.rds"
    script:
        "scripts/process_data.R"

#rule clustering
rule clustering:
    input:
        get_data_path
    output:
        clustering_plot = os.path.join(result_path, '{analysis}', 'clustering_plot_{resolution}.pdf'),
        updated_seurat_object = os.path.join(result_path, '{analysis}', 'updated_seurat_object.rds')
    params:
        resolution = config["resolution"],
    log:
        os.path.join("logs", "{analysis}_clustering.log")
    script:
        "../scripts/dea.R"

# perform differential expression analysis
rule dea:
    input:
        get_data_path
    output:
        dea_results = os.path.join(result_path,'{analysis}','DEA_results.csv'),

    resources:
        mem_mb=config.get("mem", "16000"),
    threads: config.get("threads", 1)
    conda:
        "../envs/seurat.yaml"
    log:
        os.path.join("logs","rules","dea_{analysis}.log"),
    params:
        partition=config.get("partition"),
        assay = lambda w: annot_dict["{}".format(w.analysis)]["assay"],
        metadata = lambda w: annot_dict["{}".format(w.analysis)]["metadata"],
        control = lambda w: annot_dict["{}".format(w.analysis)]["control"],
        logfc_threshold = config["logfc_threshold"],
        test_use = config["test_use"],
        min_pct = config["min_pct"],
        return_thresh = config["return_thresh"],
    script:
        "../scripts/dea.R"

# aggregate results per analysis
rule aggregate:
    input:
        dea_results = os.path.join(result_path,'{analysis}','DEA_results.csv'),
    output:
        dea_all_stats = report(os.path.join(result_path,'{analysis}','DEA_ALL_stats.csv'), 
                                  caption="../report/dea_stats.rst", 
                                  category="{}_dea_seurat".format(config["project_name"]), 
                                  subcategory="{analysis}"),
        dea_filtered_stats = report(os.path.join(result_path,'{analysis}','DEA_FILTERED_stats.csv'), 
                                  caption="../report/dea_stats.rst", 
                                  category="{}_dea_seurat".format(config["project_name"]), 
                                  subcategory="{analysis}"),
        dea_filtered_lfc = os.path.join(result_path,'{analysis}','DEA_FILTERED_LFC.csv'),
        dea_all_stats_plot = report(os.path.join(result_path,'{analysis}','plots','DEA_ALL_stats.png'), 
                                  caption="../report/dea_stats.rst", 
                                  category="{}_dea_seurat".format(config["project_name"]), 
                                  subcategory="{analysis}"),
        dea_filtered_stats_plot = report(os.path.join(result_path,'{analysis}','plots','DEA_FILTERED_stats.png'), 
                                  caption="../report/dea_stats.rst", 
                                  category="{}_dea_seurat".format(config["project_name"]), 
                                  subcategory="{analysis}"),
    resources:
        mem_mb=config.get("mem", "16000"),
    threads: config.get("threads", 1)
    conda:
        "../envs/ggplot.yaml"
    log:
        os.path.join("logs","rules","aggregate_{analysis}.log"),
    params:
        partition=config.get("partition"),
        assay = lambda w: annot_dict["{}".format(w.analysis)]["assay"],
        metadata = lambda w: annot_dict["{}".format(w.analysis)]["metadata"],
        control = lambda w: annot_dict["{}".format(w.analysis)]["control"],
        adj_pval = config["filters"]["adj_pval"],
        lfc = config["filters"]["lfc"],
        min_pct = config["filters"]["min_pct"],
        score_formula = config["score_formula"],
    script:
        "../scripts/aggregate.R"

