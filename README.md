# beetresmabs
scripts related to analyses performed for beetresmabs

# 1) VCF_combiner.py

This script can combine individual VCF files into a single VCF file suitable for the following dAF analysis.

```
Usage:
  python VCF_combiner.py --in <DIR> --out <FILE>

  --in   STR    Input folder
  --out  STR    Output VCF file
```

`--in` a folder containing VCF files that should be combined.

`--out` a VCF file that will contain one data column corresponding to each input VCF.


# 2) filter_parent_variants.py

This script can generate a gold standard VCF for the filtering of a given F2 generation VCF file. This gold standard is based on two parent VCF files and a F1 VCF file. Variants need to be homozygous in the two parents and heterozygous in the F1 generation.

```
Usage:
  python filter_parent_variants.py --p1 <FILE> --p2 <FILE> --f1 <FILE> --out <DIR>

  --p1   STR    Parent1 VCF file
  --p2   STR    Parent2 VCF file
  --f1   STR    F1 VCF file
  --out  STR    Output VCF file
```

`--p1` specifies a VCF file of the parent1 variants.

`--p2` specifies a VCF file of the parent2 variants.

`--f1` specifies a VCF file of the F1 variants.

`--out` specifies the output folder where the result files will be stored.




# 3) filter_VCF_by_goldstandard.py

# 4) dAF_selected_contigs.py

