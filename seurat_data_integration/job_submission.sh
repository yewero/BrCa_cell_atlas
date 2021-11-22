# SEURAT V3 JOB SUBMISSION SCRIPT
# this script is for the integration of seurat objects using CCA
# SUNNY WU
#
#
# Activate conda R enviroment
#source activate Renv
#
#
# RUN ID
SAMPLEIDS="CCA23Jun2019"
#
# R PATH
R="R"
TEMPPWD=$(pwd)
# output directory structure
mkdir -p output
for SAMPLENAME in ${SAMPLEIDS}; do
  mkdir -p output/CCA_${SAMPLENAME}
done
# SUBMIT JOBS
for SAMPLENAME in ${SAMPLEIDS}; do

  echo ${SAMPLENAME}
  cd output/CCA_${SAMPLENAME}

  # input sample sheet
  SUBSETINPUTFILE="${TEMPPWD}/config/sample_input_file.csv"
  # reference genome
  SPECIES="human"
  # seurat gene input file
  SEURATGENEINPUTFILE="${TEMPPWD}/config/seurat_gene_input_file.csv"
  # seurat params file
  SEURATPARAMFILE="${TEMPPWD}/config/seurat_CCA_params_file.csv"
  # seurat script
  SEURATSCRIPT="${TEMPPWD}/config/scripts/seurat_CCA_HPC_processing.R"
  # job name
  SEURATJOBNAME="sCCA_${SAMPLENAME}"
  # path to individual seurat objects
  OBJECTSPATH="${TEMPPWD}/../seurat_individual_data_processing/output"

  # cd to directory
  cd $TEMPPWD/output/CCA_${SAMPLENAME}/
  # submit job
  qsub \
  -cwd \
  -pe smp 32 \
  -l mem_requested=20G \
  -b y \
  -j y \
  -V \
  -P TumourProgression \
  -N ${SEURATJOBNAME}\
  "${R} CMD BATCH \
  --no-save '--args \
  ${SAMPLENAME} \
  ${SPECIES} \
  $(pwd) \
  ${SUBSETINPUTFILE} \
  ${SEURATGENEINPUTFILE} \
  ${SEURATPARAMFILE} \
  ${OBJECTSPATH}' \
  ${SEURATSCRIPT}"

  cd ${TEMPPWD}

  done
