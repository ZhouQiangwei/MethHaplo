## ASM analysis (BAM file with MD:Z state)
../bin/methyhaplo -M asm -a Y -m asmtestdata.mr -b asmtestdata.bam -o asmtest.output
## ASM analysis (BAM file without MD:Z state)
../bin/methyhaplo -M asm -a N -m asmtestdata.mr -b asmtestdata.nomd.bam -o asmtest.output -g hg38.chr.fa
## Hap analysis
../bin/methyhaplo -M hap -a Y -m asmtestdata.mr -b asmtestdata.bam -o asmtest.output

## results [asm bed file]
### chromosome start end length Nmeth
chr1    2189820 2189938 Meth:118        8

## results [asm detail file]
### chromosome start end MM MU UM UU pvalue qvalue
******span:     chr1    2018573 2018606 Meth:33 5
chr1    2018573 2018578 4       0       0       22      0.000067        0.010494
chr1    2018578 2018598 4       0       3       19      0.002341        0.054316
chr1    2018598 2018604 7       0       0       17      0.000003        0.002242
chr1    2018604 2018606 5       2       0       17      0.000494        0.026542

## result [hap]
### number haptype1 haptype2 chromosome pos hap1 hap2 hap1/hap2 
BLOCK: offset: 1 len: 8 phased: 8 SPAN: 511 fragments 28
1       1       0       chr1    2000047 C       T       0/1     0       .       52.60
2       1       0       chr1    2000150 C       T       0/1     0       .       100.00

### Name chromosome start end
### chromosome pos hap1 hap2 hapfile[pos/neg]
#Block  chr1    2000047 2000558
chr1    2000047 T       C       1
chr1    2000150 T       C       1
chr1    2000300 T       C       1
