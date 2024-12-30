#!/bin/bash
#SBATCH -J BM4
#SBATCH -p amd-ep2,amd-ep2-short,intel-sc3
#SBATCH -q normal
#SBATCH -c 40
#SBATCH --mem=400G

module load cellranger/6.1.1
module load R/4.2.1
module load gcc/11.2.0
module load fftw/3.3.10

### 1.preprocess
### for single-cell RNA-seq:
rawdata_path='Data/CleanData/'
sample="baseline C3D1"
for s in $sample
do
cellranger multi --id=P1_3_${s} --csv=P1_3_${s}_multi-config-template.csv --localcores=40
done
### for single-cell RNA-seq and TCR-seq:
# cellranger multi --id=D01_multi_result --csv=D01_multi-config-template.csv --localcores=40

### 2.quality control (single sample)
R_script="Rscript/"
sample="P1_3_baseline P1_3_C3D1"
rawdata_path='./'
for s in $sample
do
subpath="${s}/outs/per_sample_outs/${s}/count/sample_feature_bc_matrix/"
mkdir $s
Rscript ${R_script}/1.Seurat3_QualityControl.r -i ${rawdata_path}/1.rawdata/${subpath} -w ${rawdata_path}/2.cellanno/$s/ -o ${s} -s Human > ${rawdata_path}/2.cellanno/${s}/1.QualityControl.log
done

### 3.combined
rawdata_path='./Human/'
combined="P1_3"
mkdir $combined
cd $combined
ln -s ../../2.cellanno/P1_3_baseline/*_raw.rds P1_3_baseline_raw.rds
ln -s ../../2.cellanno/P1_3_C3D1/*_raw.rds P1_3_C3D1_raw.rds
cd ..

### PBMC_10k: downloading from https://cf.10xgenomics.com/samples/cell-vdj/5.0.0/sc5p_v2_hs_PBMC_10k/
sample_names="PBMC_10k,P1_3_baseline,P1_3_C3D1"
rds_files="PBMC_10k_raw.rds,P1_3_baseline_raw.rds,P1_3_C3D1_raw.rds"
filter_fea="6000,4000,5000"
filter_mt="8,5,5"
Rscript ${R_script}/2.IntegratedSamples.Seurat3_PCAselection.downsample10000.r -w ${rawdata_path}/3.combination/$combined -l $sample_names -f $rds_files -u $filter_fea -m $filter_mt -o ${combined} > ${rawdata_path}/3.combination/2.PCAselection.log
Rscript ${R_script}/CCA_3.Seurat3_tSNEorUMAP.r -w ${rawdata_path}/3.combination/$combined -p 30 -f ${combined}_Origin_Integrated.rds -o ${combined} > ${rawdata_path}/3.combination/3.CCA_p30.log
