# scCNAInfer
This repo provides ploidy correction given the raw read counts and copy number profiles from existing tools.
![alt text](https://github.com/compbio-mallory/SCCNAInfer/blob/main/SCCNAInfer.png)
## Contents
1. [Setup](#setup)
2. [Preprocess](#pre)
3. [Running SCCNAInfer Correction](#runSCCNAInfer)
	- [Input](#input)
	- [Output](#output)
4. [Running SCCNAInfer with Aneufinder, Ginkgo, SCOPE, SeCNV Input](#runothers)
   - [Aneufinder](#aneufinder)
   - [Ginkgo](#ginkgo)
   - [SCOPE](#scope)
   - [SeCNV](#secnv)
5. [Wrapper Script to run Aneufinder, Ginkgo, SCOPE, SeCNV with Correction](#wrapper)
6. [Miscellaneous](#mis)

<a name="setup"></a>
## Setup

The scripts are written in Python3. Following Python packages should be installed:

numpy, pandas, scipy, sklearn, linecache

The following Bioinformatic tools shoud be installed to preprocess the data. 

bwa, samtools, bedtools, bigWigAverageOverBed, pyfaidx, picard
<a name="pre"></a>
## Preprocessing
The preprocess the bam files, users can following the preprocessing steps described [here](https://github.com/deepomicslab/SeCNV).
<a name="runSCCNAInfer"></a>
## Running SCCNAInfer

<a name="input"></a>
### Input
SCCNAInfer takes outputs from an existing method and perform correction on CNV based on the clonality. Commands to run Aneufinder, Ginkgo, SCOPE, and SeCNV are provided [below](#runothers) 

When only raw read count is provided, SCCNAInfer will perform segmentation adapted from SeCNV, then perform the correction accordingly. 

`python scripts SCCNAInfer.py -cov $reads -CN $CN -path $path -gc $gc -out $out -ref $ref`
- `-cov` [Required] A tab separated read depth files. Should be in this format 'CHROM START END Cell1 Cell2 ...'
- `-path`[Required] Working path
- `-gc` [Required] GC and mappability file. Should be in this format 'CHROM START END gc mapp'. A precomuted gc and mapp files can be found in the script folder.
- `-out`[Required] Output file prefix.
- `-CN` [Optional] A tab separated segmentation file or abosolute copy number profile file. If provided, SCCNAInfer will extract the segmentation information from this file instead of performing segmentation. If segemntation file is provide, it should contains exact three columns 'CHROM START END'. If the copy number profile is provided, the file be in this format 'CHROM START END Cell1 Cell2 ...'
- `-ref` [Optional] reference type. Either hg19 or hg38. Default is hg19. 
- All input files should be in the same folder $path. 

<a name="output"></a>
### Output
SCCNAInfer will return a tsv file `$out_cnv.tsv` with corrected copy number in this format 'CHROM START END Cell1 Cell2 ...'.

<a name="runothers"></a>
## Running SCCNAInfer with Aneufinder, Ginkgo, SCOPE, SeCNV Input
<a name="aneufinder"></a>
### Run SCCNAInfer with Aneufinder Input
1. Download and install Aneufinder by following the instruction [here](https://github.com/ataudt/aneufinder)
2. Run scripts `aneufinder/run_aneufinder.R -b $bamdir -d $wd` where `$bamdir` is the directory with the bamfiles, and `$wd` is the working directory, where the output will be written into.
3. Extract the copy number profile into a formatted tsv file. `python run_aneufinder/get_aneufinderCN.py $wd/BROWSERFILES/method-edivisive/binsize_5e+05_stepsize_5e+05_CNV.bed.gz aneufinder_cnv.tsv 500000`
4. Extract the raw read counts into a formated tsv file. `
Rscript run_aneufinder/getRawReads.R binned aneufinder_reads.tsv`
5. Ready to run the correction.
6. `python scripts SCCNAInfer.py -cov aneufinder_reads.tsv -CN aneufinder_cnv.tsv -path $wd -gc gc_map.tsv -out aneufinder -ref hg19`
<a name="ginkgo"></a>
### Run SCCNAInfer with Ginkgo Input
1. Download, install and run Ginkgo by following the instruction [here](https://github.com/compbiofan/SingleCellCNABenchmark#ginkgo)
2. Copy the modified `/path/to/ginkgo/ginkgo.sh` to the `ginkgo/cli` directory. This script commented out plotting commandas.
3. Copy the modifed `/path/to/ginkgo/process.R` to the `ginkgo/scripts` directory. This script is modified in terms of output format. 
4. Get bed files. For every bam file, `bedtools bamtobed -i $basename.bam  >  base_name.bed`
5. Get the cell list. `ls *bed > cell.list` 
6. Run ginkgo, `bash /path/to/ginkgo/cli/ginkgo.sh --input $wd --genome hg19 --binning $bins --cells cell.list`
7. Format output. `python /path/to/scCNAPolish/ginkgo/addChrom2reads.py  $wd/data   $wd/SegCopy $output_Ginkgo_reads.tsv` 
8. Copy scripts/gc_ma_hg19.tsv to $wd.
9. Run correction `python /path/to/scCNAPolish/scripts/SCCNAInfer.py  -cov $out_Ginkgo_reads.tsv -CN SegCopy -path $wd -gc gc_map_hg19.tsv  -out $output  -ref $ref`

<a name="scope"></a>
### Run SCCNAInfer with SCOPE Input
MUST have normal cells otherwise the program will exit with an error. 
1. Download and install SCOPE by following the instruction [here](https://github.com/rujinwang/SCOPE/tree/master)
2. Run scripts `scope/run_scope.R -b $bamdir -d $wd -r $ref -o $out ` where `$bamdir` is the directory with the bamfiles, and `$wd` is the working directory, where the output will be written into, `ref` is the reference type either hg19, or hg38. Add `-s` if bam files are single-end. "-p $pattern" is the bam file pattern '\*.dedup.bam$'.
3. This script will return a cnv file `$out_SCOPE_cnv.tsv` and a raw read depth file `$out_SCOPE_reads.tsv`, and a gc_map.tsv file.
4. Ready to run the correction
5.  `python scripts SCCNAInfer.py -cov $out_SCOPE_reads.tsv -CN $out_SCOPE_cnv.tsv-path $wd -gc gc_map.tsv -out $out -ref $ref`
<a name="secnv"></a>
### Run SCCNAInfer with SeCNV Input
1. Get read depth file. `python scripts/preprocess.py $wd $refFile $bamdir $pattern $bin $ref`. `$wd` working directory. `$refFile` is the absolute path the the reference, eg /path/to/ref.fa.`$bamdir` is the directory holding preprocessed bam files.`$pattern` is the pattern of the preprocessed bam files, eg \*dedup.bam. `$bin` is the bin width, eg 500000. `$ref` is either hg19 or hg38. 
2. Users can use the main script to run SeCNV `python scripts/SCCNAInfer.py -cov $reads  -path $path -gc $gc -out $out -ref $ref`. SeCNV's result will be returned as intermediate results. 

To run the most recent version of SeCNV, please refer [here](https://github.com/deepomicslab/SeCNV), and cover the output CNV file into the required format 'CHROM START END Cell1 Cell2 ...' before running the correction. 

<a name="wrapper"></a>
## Wrapper Script to run Aneufinder, Ginkgo, SCOPE, SeCNV with Correction
Here we provide a wrapper script to run one of  Aneufinder, Ginkgo, SCOPE, SeCNV with correction. 

Command to run SCCNAInfer with SeCNV Input `python scCNAPolish.py -m SeCNV -bdir /path/to/dedup_bam/  -wd /path/to/wd/ -ref hg38 -fa /path/to/hg38/hg38.fa -out $output -gc gc_map.tsv -binw 500000 `

Command to run SCCNAInfer with SCOPE Input `python scCNAPolish.py -m SCOPE -bdir /path/to/dedup_bam/  -wd /path/to/wd -ref hg38 -fa /path/to/hg38/hg38.fa -out $output -gc gc_map.tsv -binw 500000 -pair`. Remove `-pair` if it is single-end. 

COmmand to run SCCNAInfer with Aneufinder Input `python scCNAPolish.py -m Aneufinder -bdir /path/to/bam_dedup/  -wd /path/to/wd/ -ref hg19 -fa /path/to/hg19/hg19.fa -out $output  -binw 500000`

Command to run SCCNAInfer with Ginkgo Input ` python scCNAPolish.py -m Ginkgo -bdir /path/to/dedup_bam/  -wd /path/to/wd/ -ref hg19 -out $output  -gbin fixed_500000 -gp /path/to/ginkgo/main/dir/`

<a name="mis"></a>
## Mis
### Calculate GC and mappability
GC and mappability are required in order to udpate reads.</br>
Use the following command to calculate GC and mappability from the segmentation from other CNV calling methods. </br>
`python cal_gc_map.py $path $ref $ref_type $bins` </br>
where `$path` is the work path. `$ref` is the path to the reference file. `$ref_type` is the reference type: hg19 or hg38. `$bins` is a tab-sep file contains (variable/fixed) consecutive bins: the first three columns should be CHROM, START, END. This will output a file `gc_map.tsv` in the `$path`. 
