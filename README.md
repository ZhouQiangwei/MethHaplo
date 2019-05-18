## Methyhaplo: Combining Allele-specific DNA Methylation and SNPs for Haplotype Assembly

## This is a README file for the usage of Methyhaplo.
------

## REQUIREMENTS
------
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
You can download the test data on ""

### Usage
        Methyhaplo: Combining Allele-specific DNA Methylation and SNPs for Haplotype Assembly\n,
        Usage: methyhaplo -M [mode] -m methfile -s <sam>/-b <bam> -o outputprefix\n,
        Options:\n,
           -M <string> [hap|asm]         methyhaolo analysis mode\n,
                                             hap: iterative approach, prefer longer haplotype results;\n,
                                             asm: hypergeometric approach, prefer accurate asm results.(default: hap);\n,
           -m, --methfile <file>         methratio file (requires)\n,
                                             format: chr  pos  strand  context  methlevel  methC  coverage\n,
           -o, --out <string>            output file prefix\n,
           -s, --sam <samfile>           sam file from batmeth2-align.  This file should be coordinate sorted, \n,
                                             using the <samtools sort> command, and must contain methylstatus[MD:Z:].\n,
           -b, --bam <bamfile>           bam file, should be coordinate sorted. (use this option or -s option but not both)\n,
           -q <int>                      only process reads with mapping quality >= INT [default >= 20].\n,
           -c, --context                 methylation context process for methyhaplo. CG, CHG, CHH, ALL[default].\n,
           -C, --NMETH                   Number of methylated reads cover cytosine site. default: 2 [m>=2]\n,
           -N, --NCOVER                  Number of coverage reads in cytosine site. default: 6 [n >= 6]\n,
	       -f, --MFloat                  Cutoff of methratio. default: 0.2 [ f =< meth <= 1-f]\n,
	       --minIS <INT>                 Minimum insert size for a paired-end read to be considered as single fragment for phasing, default 0\n,
	       --maxIS <INT>                 Maximum insert size for a paired-end read to be considered as a single fragment for phasing, default 1000\n,
	       --DBtmpsize <INT>             Maximum size of temp read store, default 12000. (only useful in asm mode)\n,
           --PE                          Paired-end reads.[default:single-end]\n,
           -v, --vcffile <file>          snp file (optional)\n,
           -r, --chromosomal-order       Use natural ordering (1,2,10,MT,X) rather then the default (1,10,2,MT,X). \n,
                                             This requires new version of the unix \sort\ command which supports the --version-sort option.\n,
           -p, --parallel <int>          Change the number of sorts run concurrently to <int>\n,
           -t, --temporary-directory     Use a directory other than /tmp as the temporary directory for sorting.\n,
           -h, -?, --help                This help message.\n,

