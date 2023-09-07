# scPloidy

### Calculate GC and mappability
GC and mappability are required in order to udpate reads.</br>
Use the following command to calculate GC and mappability from the segmentation from other CNV calling methods. </br>
`python cal_gc_map.py $path $ref $ref_type $bins` </br>
where `$path` is the work path. `$ref` is the path to the reference file. `$ref_type` is the reference type: hg19 or hg38. `$bins` is a tab-sep file contains (variable/fixed) consecutive bins: the first three columns should be CHROM, START, END. This will output a file `gc_map.tsv` in the `$path`. 
