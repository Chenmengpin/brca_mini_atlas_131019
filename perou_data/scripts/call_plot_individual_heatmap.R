#!/bin/bash

#sample_names=( HTP196 HTP196_epithelial_control WHIM11 HT35A4AL HT3821AL CID44971)
sample_names=( HTP196 WHIM11 )

project_name="identify_epithelial"
subproject_name="brca_mini_atlas_131019"
ncores=10
cancer_x_threshold_sd_multiplier=2
cancer_y_threshold_sd_multiplier=1.5
normal_x_threshold_sd_multiplier=1
normal_y_threshold_sd_multiplier=1.25
include_group_annotation="FALSE"
adjust_normal_thresholds="FALSE"

home_dir="/share/ScratchGeneral/jamtor"
project_dir="$home_dir/projects/single_cell/$project_name/$subproject_name/perou_data"
script_dir="$project_dir/scripts/"

for sample_name in ${sample_names[@]}

  do echo "Running InferCNV for $sample_name..."

  log_dir="$project_dir/logs/$sample_name"
  mkdir -p $log_dir
  echo "Logs are in $log_dir"
  qsub -wd $log_dir -pe smp $ncores -N icnv.$sample_name -b y -j y -V -P TumourProgression \
    "/share/ClusterShare/software/contrib/CTP_single_cell/tools/R_developers/config_R-3.5.0/bin/R \
    CMD BATCH  --no-save '--args $sample_name $ncores $cancer_x_threshold_sd_multiplier $cancer_y_threshold_sd_multiplier $normal_x_threshold_sd_multiplier $normal_y_threshold_sd_multiplier $include_group_annotation $adjust_normal_thresholds' $script_dir/plot_individual_heatmap.R"
  #  /share/ClusterShare/software/contrib/CTP_single_cell/tools/R_developers/config_R-3.5.0/bin/Rscript --vanilla $script_dir/plot_individual_heatmap.R  $sample_name $ncores $cancer_x_threshold_sd_multiplier $cancer_y_threshold_sd_multiplier $normal_x_threshold_sd_multiplier $normal_y_threshold_sd_multiplier $include_group_annotation $adjust_normal_thresholds

  printf "\n"
done
