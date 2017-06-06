# biodatapooler

## Table of contents
1. Python version
2. Dependencies
3. Usage
4. Output
5. Post-processing

### 1. Python version 
Python 2.7

### 2. Dependencies
Python modules: numpy1.7.1, BioPython (AlignIO, Align, Phylo), scipy,
translatorx (http://translatorx.co.uk/),
PAML (http://abacus.gene.ucl.ac.uk/software/),
Perl

### 3. Usage

General usage:
```
python -u bioadatapooler.py -c controlfile.txt >  output.txt
```
This is a data integration pipeline with a lot of subroutines that are designed to 
integrate data for sets of ortholgous proteins accross multiple species. 
Input data is presumed to be in form of files of various formats (fasta, gff3, cuffdiff, 
for now), each containing information for only one species. 
All data files are to be listed in the control file that is supplied to this code. 
User also supplies the address of the working directory and control files for running PAML 
using the control file. In future, paths of all external tools being used here (PAML, 
translatorx, maftools etc.) will be supplied using this control file as well.
Control file has a tab-separated format that needs to be followed strictly. 
Please see the supplied example files for proper usage.

Example:

```
#speciesTag	type	method	file
ALL	numspecies	param	4
ALL	orthomcl	ortholist	orthoMCL.output.file.txt
ALL	pamlctl	pairwise	pairwise_paml_controlfile.ctl
ALL	pamlctl	null	paml_ModelM0_controlfile.ctl
ALL	pamlctl	freeRatio	paml_FreeRatio_controlfile.ctl
MUNI	pamlctl	neutral	paml_MuniBranchNeutral_controlfile.ctl
MUNI	pamlctl	special	paml_MuniBranchSpecial_controlfile.ctl
ALL	workdir	NA	your/working/directory/
MRAP	annotation	gff	MRAP.gff3
MUNI	annotation	gff	MUNI.gff3
MUNI	cuffdiff	male_female	gene_exp.diff
MUNI	annotation	protein	Muni.proteins.fasta
MUNI	annotation	nucleotide	Muni.nucleotides.fasta
MUNI	annotation	genome	muni.genome.fasta
```

For the control file, valid 'type' and related 'method' are:

```
TYPE          METHOD
numSpecies    param
orthoml       ortholist
pamlctl       pairwise,null,freeRatio,neutral,special
workdir       NA
annotation    gff,protein,nucleotide,genome
columndata    NA
```

This list will expand, allowing users to enter paths of external tools. 

Integration of these
species-specific datasets is driven by a user-provided list of proteins; most likely
this list is of orthologous proteins. The ortholgous protein list is of the format of
OrthoMCL output:
<name of orthologous group>: <species tag>|<protein 1> <species tag>|<protein 2> ...

For example
MN123: hg19|abc1 hg19|abc2 mm|abc1 mm|abc3 
The orthology file is listed in control file as well.
This orthology information is loaded into a special class, which has multiple
subroutines of itself, and which connects to other subroutines that together can perform
multiple tasks for the user:
1. align at nucleotide level
2. align at amino acid sequence level
3. obtain pair-wise dN/dS estimates using PAML, from multiple sequence alignment
4. obtain lineage-specific dN/dS estimates
5. perform log-likelihood tests of selection and neutral evolution using PAML


###4. Output
The output in a 'melted' data format, in the jargon of R's reshape2 file formats. Every data point accumulated is listed in a separarate line, with all the associated information for that data poined listed in the preceeding tab-separated columns. 

Roughly, the columns of the following example correspond to this:

unused   orthologous    tool   condition valueType     value
column    group name     used   
```
#data   group         analysis  species  variable        value
data    MB_PB00083_C    paml    null    start:kappa     2.93277
data    MB_PB00083_C    paml    null    start:dN        0.0186
data    MB_PB00083_C    paml    null    start:codon model       F3x4
data    MB_PB00083_C    paml    null    start:numParams 9
data    MB_PB00083_C    paml    null    start:dS        0.5903
data    MB_PB00083_C    paml    null    start:codon__model      F3x4
data    MB_PB00083_C    paml    null    start:version   4.8a
data    MB_PB00083_C    paml    null    start:omega     0.03157
data    MB_PB00083_C    paml    null    start:lnL       -2342.386929
data    MB_PB00083_C    paml    null    start:tree length       0.45407
data    MB_PB00083_C    paml    null    start:description       one-ratio
data    MB_PB00083_C    paml    null    start:model     One dN/dS ratio for branches,
data    MB_PB00083_C    paml    MUNI.special    7..1:dN 0.0007
data    MB_PB00083_C    paml    MUNI.special    start:dN        0.0186
data    MB_PB00083_C    paml    MUNI.special    6..7:S  312.6
data    MB_PB00083_C    paml    MUNI.special    5..6:N*dN       8.7
data    MB_PB00083_C    paml    MUNI.special    6..7:N  1034.4
data    MB_PB00083_C    paml    MUNI.special    5..4:S  312.6
data    MB_PB00083_C    paml    MUNI.special    6..7:S*dS       8.6
```

###5. Post-processing
 
The user is expected to use R's reshape2 package to summarize this information in form of a matrix, that could be used for plotting and statistical analysis. The reshape command used for above output file is:

```
d<-read.table("<input file>", header=T, sep="\t") 
library(reshape2) 
aqw <- dcast(d, group ~ analysis + species + variable , value.var="value") 
write.table(aqw, "Reshaped_data.txt", quote=F, sep="\t", row.names=F) 
```




