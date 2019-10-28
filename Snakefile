# snakemake --reason --cores 7 --cluster "qsub -pe smp 7 -N smk.infercnv -b y -j y -V -P TumourProgression" -j 23

configfile: "config.yaml"

# directories:
R_dir = "/share/ClusterShare/software/contrib/CTP_single_cell/tools/R_developers/config_R-3.5.0/bin"

rule all:
    input:
        expand(
            "results/infercnv/t_cells_included/{sample}/" +
            "threshold_sd_multipliers_{cancer_x_threshold_sd_multiplier}_" +
            "{cancer_y_threshold_sd_multiplier}_" +
            "{normal_x_threshold_sd_multiplier}_" +
            "{normal_y_threshold_sd_multiplier}/" +
            "plots/CNA_density_plot.png",
            sample=config["samples"],
            cancer_x_threshold_sd_multiplier=config["cancer_x_threshold_sd_multiplier"],
            cancer_y_threshold_sd_multiplier=config["cancer_y_threshold_sd_multiplier"],
            normal_x_threshold_sd_multiplier=config["normal_x_threshold_sd_multiplier"],
            normal_y_threshold_sd_multiplier=config["normal_y_threshold_sd_multiplier"]
        )

rule infercnv:
    input:
        "raw_files/seurat_objects/{sample}/" +
            "03_seurat_object_processed.Rdata"
    output:
        final_observations = protected(
            "results/infercnv/t_cells_included/{sample}/" +
            "infercnv.12_denoised.observations.txt"
        )
    params:
        subset_data=config["subset_data"],
        include_t_cells=config["include_t_cells"]

    threads: 6
    shell:
            "mkdir -p logs/{wildcards.sample}/; " +
            "cd logs/{wildcards.sample}/; " +
            R_dir + "/R CMD BATCH  --no-save '--args {wildcards.sample} " +
            "{params.subset_data} {params.include_t_cells}' ../../scripts/1.infercnv.R"

rule individual_plot:
    input:
        "results/infercnv/t_cells_included/{sample}/" +
        "infercnv.12_denoised.observations.txt"
    output:
        infercnv_plot = "results/infercnv/t_cells_included/{sample}/" +
            "threshold_sd_multipliers_{cancer_x_threshold_sd_multiplier}_" +
            "{cancer_y_threshold_sd_multiplier}_" +
            "{normal_x_threshold_sd_multiplier}_" +
            "{normal_y_threshold_sd_multiplier}/" +
            "plots/CNA_density_plot.png",
        density_plot = "results/infercnv/t_cells_included/{sample}/" +
            "threshold_sd_multipliers_{cancer_x_threshold_sd_multiplier}_" +
            "{cancer_y_threshold_sd_multiplier}_" +
            "{normal_x_threshold_sd_multiplier}_" +
            "{normal_y_threshold_sd_multiplier}/" +
            "plots/infercnv_plot.png",
        quad_plot = "results/infercnv/t_cells_included/{sample}/" +
            "threshold_sd_multipliers_{cancer_x_threshold_sd_multiplier}_" +
            "{cancer_y_threshold_sd_multiplier}_" +
            "{normal_x_threshold_sd_multiplier}_" +
            "{normal_y_threshold_sd_multiplier}/" +
            "plots/normal_call_quad_plot_mean_of_scaled_squares.png",
    params:
        include_t_cells=config["include_t_cells"],
        cancer_x_threshold_sd_multiplier=config["cancer_x_threshold_sd_multiplier"],
        cancer_y_threshold_sd_multiplier=config["cancer_y_threshold_sd_multiplier"],
        normal_x_threshold_sd_multiplier=config["normal_x_threshold_sd_multiplier"],
        normal_y_threshold_sd_multiplier=config["normal_y_threshold_sd_multiplier"]
    threads: 6
    shell:
        "cd logs/{wildcards.sample}/; " +
        R_dir + "/R CMD BATCH  --no-save '--args {wildcards.sample} " +
        "{params.include_t_cells}  " +
        "{params.cancer_x_threshold_sd_multiplier} " +
        "{params.cancer_y_threshold_sd_multiplier} " +
        "{params.normal_x_threshold_sd_multiplier} " +
        "{params.normal_y_threshold_sd_multiplier}' ../../scripts/2.plot_individual_heatmap.R"
