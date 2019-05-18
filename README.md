## Methyhaplo: Combining Allele-specific DNA Methylation and SNPs for Haplotype Assembly

DNA methylation is an important epigenetic modification that plays a critical role in most eukaryotic organisms. Parental alleles in haploids may exhibit different methylation patterns, which can lead to different phenotypes, and even different therapeutic and drug responses to diseases. However, there is currently no suitable software to obtain accurate DNA methylation haplotype results to our knowledge. To address this issue, we developed a new method, Methyhaplo, for haplotype assembly with allele-specific DNA methylation and SNPs from whole-genome bisulfite sequencing (WGBS) data. Our results showed that the haplotype assembly with allele-specific DNA methylation and SNPs was ten times longer than that with just SNPs. Moreover, Methyhaplo could integrate WGBS-Seq and Hi-C to obtain the better haplotype results. 

## This is a README file for the usage of Methyhaplo.
------

## REQUIREMENTS
1. gcc (v4.8) , gsl library
2. SAMtools
3. Python3
4. Perl

## INSTALL
------
a) Download
`git clone https://github.com/ZhouQiangwei/Methyhaplo.git`

b) Change directory into the top directory of Methyhaplo
`cd Methyhaplo`

c) Type
- make
- make install

d) The binary of Methyhaplo will be created in current folder

## USAGE of Methyhaplo
------
### Example data
You can download the test data on 

### Usage
1. MethyHaplo command
```
        Methyhaplo: Combining Allele-specific DNA Methylation and SNPs for Haplotype Assembly
        Usage: methyhaplo -M [mode] -a Y -m methfile -s <sam>/-b <bam> -o outputprefix
        Options:
           -M <string> [hap|asm]         methyhaolo analysis mode
                                             hap: iterative approach, prefer longer haplotype results;
                                             asm: hypergeometric approach, prefer accurate asm results.(default: hap);
           -m, --methfile <file>         methratio file (requires)
                                             format: chr  pos  strand  context  methlevel  methC  coverage
           -o, --out <string>            output file prefix
           -s, --sam <samfile>           sam file from batmeth2-align.  This file should be coordinate sorted, 
                                             using the <samtools sort> command, and must contain methylstatus[MD:Z:].
           -b, --bam <bamfile>           bam file, should be coordinate sorted. (use this option or -s option but not both)
           -a <Y/N>                      If bam/sam file contain MD state by batmeth2 calmeth scripts.
                                             If not, please define genome location by -g paramater.
           -g, --genome <genome>         If bam/sam file isnot contain MD.
           -q <int>                      only process reads with mapping quality >= INT [default >= 20].
           -c, --context                 methylation context process for methyhaplo. CG, CHG, CHH, ALL[default].
           -C, --NMETH                   Number of methylated reads cover cytosine site. default: 2 [m>=2]
           -N, --NCOVER                  Number of coverage reads in cytosine site. default: 6 [n >= 6]
	       -f, --MFloat                  Cutoff of methratio. default: 0.2 [ f =< meth <= 1-f]
	       --minIS <INT>                 Minimum insert size for a paired-end read to be considered as single fragment for phasing, default 0
	       --maxIS <INT>                 Maximum insert size for a paired-end read to be considered as a single fragment for phasing, default 1000
	       --DBtmpsize <INT>             Maximum size of temp read store, default 12000. (only useful in asm mode)
           --PE                          Paired-end reads.[default:single-end]
           -v, --vcffile <file>          snp file (optional)
           -r, --chromosomal-order       Use natural ordering (1,2,10,MT,X) rather then the default (1,10,2,MT,X). 
                                             This requires new version of the unix \sort\ command which supports the --version-sort option.
           -p, --parallel <int>          Change the number of sorts run concurrently to <int>
           -t, --temporary-directory     Use a directory other than /tmp as the temporary directory for sorting.
           -h, -?, --help                This help message.
```

2. Allele-specific DNA methylation region visulization
```bash
python methpoint.py chr18_26621440_26621650.md.sort.bam chr18:26621570-26621650 . chr18_26621570_26621650.UM 0
```

<p align="center">
        <img src="scripts/asmexample.png" alt="asmexample"  width="600" height="500">
</p>



