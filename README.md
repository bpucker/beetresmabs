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
Result files: P1_hom.vcf contains the homozygous variants of parent1, P2_hom.vcf contains the homozygous variants of parent2, F1_het.vcf contains the heterozygous variants in the F1, and gold_standard.vcf contains the variants supported by this analysis.

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

This script filters a given VCF file based on the goldstandard VCF. Only variants at a position listed in the goldstandard are kept.

```
Usage:
  python filter_VCF_by_goldstandard.py --in <FILE> --gold <FILE> --out <DIR>

  --in    STR    Input VCF file
  --gold  STR    Goldstandard VCF file
  --out   STR    Output VCF file
```

`--in` specifies the input VCF file.

`--gold` specifies a goldstandard VCF file that is used for the filtering.

`--out` specifies the output VCF file.


# 4) dAF_selected_contigs.py

This script performs a delta allele frequency analysis.

```
Usage:
  python dAF_selected_contigs.py --input_vcf <FILE> --reference_file <FILE> --output_dir <DIR> --pool1 <STR> --pool2 <STR>

  --in     STR    Input VCF file
  --ref    STR    Reference FASTA file
  --out    STR    Output directory
  --pool1  STR    Samples of pool1
  --pool2  STR    Samples of pool2
```

`--input_vcf` specifies the input VCF file.

`--reference_file` specifies a FASTA file containing the reference genome sequence.

`--output_dir` specifies the output directory.

`--pool1` comma-separated list of samples that form pool1.

`--pool2` comma-separated list of samples that form pool2.


# References

This repository.
