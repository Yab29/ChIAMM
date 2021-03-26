## **ChIAMM**
## Introduction
**Ch**romatin **I**nteraction **A**nalysis with **P**aired-**E**nd **T**ag (**ChIA-PET**) sequencing is a technology to study genome-wide long-range chromatin interactions bound by protein factors. **ChIAMM** is a statistical technique design for detecting significant chromatin interactions from the ChIA-PET experiments data using the Mixture model. 
### The ChIAMM requires the following dependencies:
> 1) &nbsp; R (≥3.4.0) <br />
> 2) &nbsp; StanHeaders (≥ 2.18.1)<br />
> 3) &nbsp; rstan (version ≥ 2.18.2)<br />
> 4) &nbsp; ggplot2 (≥ 2.0.0) <br />
> 5) &nbsp; tidyverse (≥ 1.2.1)<br />
> 6) &nbsp; bayesplot ( ≥1.7.1)<br />

### Data preparation
Before executing the **ChIAMM**, you need to analyze the raw ChIA-PET data using ChIA-PET Tool V3 without any FDR cutoff value. From the ChIA-PET Tool output files, the file out.cluster.withpvalue.txt will be used for downstream analysis in ChIAMM. Specifically, the first 11 columns in the file out.cluster.withpvalue.txt will be used: chrom1, start1, end1, chrom2, start2, end2, ipet count, type, distance, tag count within anchor 1 and tag count within anchor 2. We call chrom1, start1, end1 as Anchor1, and chrom2, start2, end2 as Anchor2.  
### Computing the Systematic Biases 
####  1) Average tag counts 
Call the variables tag count within anchor 1 and tag count within anchor 2 as tagcou1 and tagcou2 respectively, and compute the average of it call the variable name “tagcouAvg”. 
#### 2) Self-ligation PETS
Using the out.spet file in the ChIA-PET Tool V3 output, we can compute the self-ligation PETs of anchors using the commands as follows:  
> 1) &nbsp; awk '{if($2<$5){print $1"\t"$2"\t"$5}else{print $1"\t"$5"\t"$2}}' out.spet  > out.spet.bed3 <br />
> 2) &nbsp; bedtools coverage -a Anchor1.bed -b out.spet.bed3| cut -f4 > self1.bed <br />
> 2) &nbsp; bedtools coverage -a Anchor2.bed -b out.spet.bed3| cut -f4 > self2.bed <br />
> 2) &nbsp; then compute the average of it and call the variable name “selfAvg”  <br />

#### 3) Mappability 
Download the mappability score for human genome version hg19 using http://genome-asia.ucsc.edu/cgi-bin/hgFileUi?db=hg19&g=wgEncodeMapability. The file is in bigWig format, and we need to convert it to .bed format using bigWigToWig and BEDOPS wig2bed. Then, find the overlap region between .bed file and Anchor regions using Bedtools map. 
> 1) &nbsp; bedtools map -a Anchor1.bed -b mappability.bed -c 5 -o mean| cut -f4 > mappability1.bed <br />
> 2) &nbsp; bedtools map -a Anchor2.bed -b mappability.bed -c 5 -o mean| cut -f4 > mappability2.bed <br />
> 3) &nbsp; Compute the average of mappability1.bed and mappability2.bed, and call the variable name “mappaAvg” <br />
> 4) &nbsp; If you don’t find the mappability in the above link, you can prepare by yourself using ngs-tools.  <br />
- Example for rice RS1 reference genome we can find like
  - sh mappability.sh -i mhRS63.fa -l 35 -p mappability > mappability.bed 

#### 4) GC content 
The GC content of anchors are computed using bedtools nuc <br />
> 1) &nbsp;	bedtools nuc -fi  hg19.fa -bed  Anchor1.bed| cut -f5  > gc1.bed <br />
> 2) &nbsp;	bedtools nuc -fi  hg19.fa -bed  Anchor2.bed| cut -f5  > gc2.bed <br />
> 3) &nbsp; Compute the average of gc1.bed and gc2.bed and call variable name “gcAv” <br />

### Input files
We need to organize the input file for inter- and intra-chromosomal interaction data separately. As an example, the input file of intra-chromosomal interactions should have such kind of variables arrangement.   

   |chrom1 |start1|end1  |chrom2 |start2 |end2  |ipet  |distance|tagcou1 |tagcou2  |tagcouAvg |selfAvg|gcAv  |mappaAvg|
   |-------|-------|------|-----|-------|------|------|--------|--------|---------|----------|-------|------|--------|
   |chr1	  |840068 |840732|chr1 |855579|	856565|	2    |15672   |	3|	6|	4.5|	9.0|	0.66|	0.69|
   
### Test data sets
> 1) &nbsp;	RNAPII ChIA-PET data from human MCF7:<br /> https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM832458 and https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM832459
> 2) &nbsp;	RNAPII ChIA-PET data from human K562:<br /> https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM832464 and https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM832465 
> 3) &nbsp; CTCF ChIA-PET data from human MCF7:<br /> https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM970215 
> 4) &nbsp;	CTCF ChIA-PET data from human K562:<br /> https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM970216 
> 5) &nbsp;	RNAPII ChIA-PET data from rice MH63: <br /> https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM3767548 and https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM3767549 
> 6) &nbsp;	H3K9me2 ChIA-PET data from rice MH63:<br /> https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM3767546 and https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM3767547 

### Usage
	Rscript ChIAMM.R -h
  	Usage: ChIAMM.R [-[-input|i] <character>] [-[-prefix|p] [<character>]] [-[-inter|e]] [-[-iter|r] [<integer>]] [-[-warmup|w] [<integer>]] [-[-help|h]]
         -i| --input      input file
         -p| --prefix     output prefix (default "out")
         -e| --inter      the input is inter-chromosomal, if not, intra-chromosomal interactions data
         -r| --iter       iteration number (default 5000)
         -w| --warmup     warmup value (default 750)
         -h |--help       print help
### Example 
For intra- and inter-chromosomal interactions, we need to run the ChIAMM.R separately.
-	For intra-chromosomal interaction data:
     - ChIAMM.R -i intra_input_file.txt 
- For inter-chromosomal interaction data 
     - ChIAMM.R -i inter_input_file.txt -e    
 
### Result file
We will get the result file names out_significant_interaction.txt. 

|chrom1|start1|	end1|	chrom2|	start2|	end2|	ipet|	W1i|
|------|------|------|-----|-------|------|------|-----|
|chr1|	919077|	919880|	chr1|	998944|	999920|	4|	0.59|

- **chrom1:** The name of the chromosome on which the cluster anchor 1 exists
- **start1:** The start coordinates of cluster anchor 1
- **end1:** The end coordinate of cluster anchor 1
- **chrom2:** The name of the chromosome on which the cluster anchor 2 exists
- **start2:** The start coordinates of cluster anchor 2
- **end2:** The end coordinate of cluster anchor 2
- **ipet:** Number of PETs between cluster anchor 1 and cluster anchor 2
- **W1i :** The probability of pair 𝑖 being a true pair

##### &nbsp;CopyRight &#169; 2021 [Guoliang's Lab](http://guolianglab.org/subpages/REASERACH/publications.php), All Rights Reserved

