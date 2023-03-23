# EagleImp
## _A combined genotype phasing and imputation tool (with optional FPGA acceleration)_
EagleImp combines genotype phasing and imputation in a single tool. The algorithms are based on 'Eagle2' by Po-Ru Loh et al. _(Nature Genetics, 2016, https://pubmed.ncbi.nlm.nih.gov/27694958/)_ and 'PBWT' by Richard Durbin _(Bioinformatics, 2014, https://pubmed.ncbi.nlm.nih.gov/24413527/)_. For EagleImp, we have put a lot of effort in modifying the methods and implementation to achieve higher efficiency and faster processing. This includes the use of multi-processing features (especially in the imputation process), which makes the tool well suited for workstations or high-performance computers with many cores and lots of RAM, but common desktop computers or laptops also benefit from faster execution times.

The EagleImp paper is published in *Bioinformatics* and available online at 
https://academic.oup.com/bioinformatics/article/38/22/4999/6706779
. 

EagleImp also adds support for FPGA-based accelerated phasing if you have an Alpha Data ADM-PCIE-8K5 FPGA accelerator card available. The release of the FPGA part is in preparation. Please contact us if you want to use EagleImp with FPGA acceleration.

### Prerequisites
EagleImp has been tested on an Ubuntu 22.04 Linux system, but should also work on similar (especially later) distributions.
Compilation depends on the *development* files of several system libraries, in particular **OpenMP**, **zlib**, **BOOST** (filesystem and program_options), **TBB2** and others.
On Ubuntu, most dependencies can already be resolved by the following installation:
```
$ sudo apt install zlib1g-dev libboost-dev libboost-filesystem-dev libboost-program-options-dev libtbb2-dev
```
It also depends on the installation of **HTSlib** to read and write VCF files (<https://github.com/samtools/htslib>). Please follow the instructions listed in the file **INSTALL** there to install **HTSlib**.

If you want to take advantage of our automation script `launch_eagleimp` for simultaneous analysis of multiple input files, including automated chromosome X splitting into pseudo-autosomal regions (PARs), you need to have *GNU Awk* (**gawk**, included in most Linux distributions) or similar and **bcftools** installed as well (<https://github.com/samtools/bcftools>).
You also need the basic tools for compiling C++ code, i.e. the **build-essential** package and a recent *GCC* (**g++**).
Since EagleImp is a **CMake** project, you must have the **cmake** package installed:
```
$ sudo apt install cmake
```

### Building the executable
Assuming you want to download (or clone) the source code to a local directory "software" in your home folder:
```
$ cd ~/software
$ git clone git@github.com:ikmb/eagleimp.git
```
Change into the source code folder and build the make files:
```
$ cd eagleimp
$ cmake -DCMAKE_BUILD_TYPE=Release -S src -B build
```
This should generate a "build" folder for compiling the "Release" version of EagleImp.
If you get an error message, you probably forgot to install a dependency package.
After successfully generating the build files, change into the "build" folder and type "make" (or "make -j" on multi-core systems):
```
$ cd build
$ make -j
```
After successful compilation, the `eagleimp` executable file is now in this folder. You can (optionally) link it to a location included in the "PATH" variable of your shell, e.g. "~/bin":
```
$ ln -s eagleimp ~/bin/eagleimp
```
Try to launch. (The following only works if you linked the executable to your "PATH", otherwise call it directly, e.g. "./eagleimp"):
```
$ eagleimp -v
--------------------------------------------------------
--- EagleImp: Genotype Phasing + Imputation
---
--- v1.00_main_ed6bc8a
--- (compiled on Nov 15 2021, 16:51:35)
---
--- Copyright (C) 2018-2021 by Lars Wienbrandt,
--- Institute of Clinical Molecular Biology, Kiel University
--- Distributed under the GNU GPLv3 license.
--------------------------------------------------------

Provided arguments:
  -v
This is version v1.00_main_ed6bc8a, compiled on Nov 15 2021, 16:51:35
Send bugs to Lars Wienbrandt <l.wienbrandt@ikmb.uni-kiel.de>
```

### Quick start
We have created a minimal example that you can find in the _example_ folder in this repository. To run the example, simply type:
```
$ cd example
$ eagleimp --geneticMap genetic_map_hg19_chr22.txt --ref 22.example.ref.vcf.gz --target 22.example.vcf.gz
```
However, the provided data is only exemplary. For a real analysis you need the following:
1. To use EagleImp you need a **genetic map**. For human genome builds GRCh37 and GRCh38 they can be found e.g. at our server
<https://hybridcomputing.ikmb.uni-kiel.de/downloads>.
The map is split by chromosomes which is recommended for a faster runtime, but a single file for the whole genome would also do.<br>
**NOTE:** If your genetic map file contains a chromosome identifier column, this column **must not** contain any literals, such as "chr" or "X". Please use the corresponding numerical value instead, e.g. "23".

2. Phasing and imputation are always performed against a **haplotype reference**. A prominent and freely available reference is the **1000 Genomes** reference panel at
<ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/>.<br>
Reference files can be either in _.vcf.gz_, _.bcf_ or in our own _.qref_ format. A reference file must not contain data from more than one chromosome. So, ideally, you should have one reference file ready per chromosome.
If you are using a _.vcf.gz_ or _.bcf_ reference file, the file must be **indexed** using either _tabix_ or _bcftools_.<br>
**NOTE:** EagleImp requires the reference file be named after the chromosome it contains at the beginning of its filename, i.e. the **filename must start with the chromosome number** (or chromosome name in case of X and Y, although 23 and 24 is also allowed). A preceding "chr" is also allowed.<br>
**NOTE:** If you use a VCF/BCF reference file, the chromosome encodings in your target and the reference must be the same. This is a restriction imposed by the HTSlib reader (that also applies to the original Eagle2). We recommend using a Qref that eliminates this restriction (see _Qref format_ below). We have prepared the 1000 Genomes GRCh37/hg19 release in _.qref_ format available for download at <https://hybridcomputing.ikmb.uni-kiel.de/downloads>.

3. Finally, you need your target data in an **indexed** _.vcf(.gz)_ or _.bcf_ file. 
Again, a target file must not contain data from more than one chromosome. Please note the special instructions for chromosomes X and Y below if you want to impute data from these chromosomes.

Phasing and imputation are then simply started as in the example:
```
$ eagleimp --geneticMap <path/to/genetic_map.file> --ref <path/to/chrN.reference.file.vcf.gz/.qref> --target <path/to/target.vcf.gz>
```
After successful imputation you will find a file <path/to/target.imputed.vcf.gz> in your target folder. The *.vcf.gz* format is the default output format. A different format, e.g. *.bcf*, can be specified with the "--vcfOutFormat" option, and a different output prefix with the "-o" option. Please see "eagleimp --help" for more information about the available options.

**NOTE:** To index your file with *bcftools* simply type:
```
$ bcftools index <path/to/your.vcf.gz>
```

### Whole genome analysis
In order to simplify the process of a whole genome analysis, we provide a bash script `launch_eagleimp` that automatically starts the imputation for all files in a given folder. The script also allows you to take advantage of your computing system's multi-processing capabilities, as it can perform the analysis of multiple files in parallel.

The script also takes care of correct partitioning of chromosome X into PAR/nonPAR regions and indexes files not previously indexed before analysis. You must have *bcftools* and *GNU Awk* (*gawk*) installed for the script to work correctly.
Please have a Qref reference file (see "Qref format" below) and a genetic map ready for each target chromosome you wish to analyse. This includes reference files for the specific PAR/nonPAR regions of chromosome X if this is one of your targets. In this case, please see the notes on _Phasing and Imputation of chromosome X_ below.

**NOTE:** The filenames of your targets **have to** start with the chromosome number (or identifier in the case of X or Y) followed by a dot ".", as the script automatically chooses the correct reference file and genetic map according to the filename of a target. A preliminary "chr" is allowed as well as the following PAR-identifiers for chromosome 23/X: 23_nonPAR, 23_PAR1, 23_PAR2, X_nonPAR, X_PAR1, X_PAR2. For example, allowed target filenames could be "1.vcf.gz", "chr22.noATCG.vcf.gz", "chr23_nonPAR.bcf", "X_PAR1.mystudy.bcf", "chrX.vcf.gz"

1. `launch_eagleimp` is located in the main folder of this repository. It **must be configured** before first use. We recommend copying the script to a location that is included in the "PATH" of your shell, e.g. your binary home folder *~/bin/*. Please then follow the instructions in the script file to configure it.

2. If the script is configured correctly, you can run it with the folder containing your target files as the first argument. Additional arguments are provided transparently to the EagleImp executable.<br>
**NOTE:** The script automatically sets the arguments for target, reference and genetic map files, along with the options for multi-threading and main memory partitioning. It also adds the default options that you configure in the script.

3. In addition to the regular output of each EagleImp run, the script generates a log file for each run and creates a summary of all runs.

##### Example
To run EagleImp on your target VCF's and add the option "\-\-outputPhasedFile" to the runs in order to produce an additional phased file of your input targets, simply type:
```
$ launch_eagleimp <your_target_folder> --outputPhasedFile
```

### Qref format
Since loading a reference panel from a VCF format can take a long time, we recommend creating a Qref file from a VCF reference if phasing and imputation is to be performed more than once with the same reference panel.
Using a Qref rather than a VCF reference has more advantages: HTSlib, which is also used in original Eagle2 to read VCF files, does not recognize two chromosomes as the same in target and reference if one is encoded only with the chromosome number and the other prefixed with "chr", e.g. _22_ and _chr22_ are not the same for the HTSlib reader. In contrast, it does not matter if you created a Qref reference from a VCF file where the chromosome is encoded with or without leading "chr", and it does not matter which encoding is present in your target file. The same applies to the encoding of chromosomes X and Y which can be encoded either literally with the character "X" or "Y" or as the number "23" or "24", respectively. In addition, the prefix "chr" may or may not be present.

Creating a ".qref" is simple:
```
eagleimp --ref <path/to/chrN.reference.file.vcf.gz> --makeQref
```
The created ".qref" file can now be used for phasing and imputation by replacing the VCF reference at the "\-\-ref".
For the minimal example it looks like this:
```
$ cd example
$ eagleimp --ref 22.example.ref.vcf.gz --makeQref
$ eagleimp --geneticMap genetic_map_hg19_chr22.txt --ref 22.example.ref.qref --target 22.example.vcf.gz
```

**NOTE:** The option "\-\-noMissingIDs" is automatically enabled when creating a Qref. This converts missing reference IDs to a _chr:position:ref:alt_ format.

**NOTE:** As already mentioned above, we have prepared the 1000 Genomes GRCh37/hg19 release in _.qref_ format available for download at <https://hybridcomputing.ikmb.uni-kiel.de/downloads>.

### The most important options
The following is a brief introduction to some options that may be useful in many use cases. For a complete list of options, see:
```
$ eagleimp --help
```
##### \-\-excludeMultiAllRef
This option excludes all variants that are multi-allelic in the reference. Otherwise, multi-allelic variants will be split into multiple bi-allelic variants where haplotypes will be set to the reference-allele if the allele is not included in the new bi-allelic variant.

##### \-\-allowRefAltSwap
It often happens that the alleles of the target variants are flipped with respect to the reference, i.e. the reference allele of a variant in the target is the alternative allele in the corresponding reference variant. By default, such variants are not considered to be shared and handled separately: In this case, the target variant is discarded as "not found in the reference", but imputed from the reference (as a reference-only variant), with the notation of the alleles as in the reference. With the "\-\-allowRefAltSwap" switch, such variants are considered shared between target and reference, and the variant will be phased and is used as an anchor in the imputation.
For example, if your target variant has alleles A and G, but the corresponding reference has G and A at the same position, your target variant will also be converted to G and A, and all haplotype codes will be inverted.

**NOTE:** The output contains the alleles as they appear in the reference.

##### \-\-allowStrandFlip
Similar to the situation above, it sometimes happens that variants are genotyped from the opposite strand. Such variants are considered shared when the "\-\-allowStrandFlip" switch is enabled.
For example, if your target variant has alleles A and G, but the corresponding reference has T and C at the same position, your target variant will also be converted to T and C, as this is considered the same variant typed from the opposite strand.

**NOTE:** The output will contain the alleles as they appear in the reference.

**NOTE:** If used in combination with "\-\-allowRefAltSwap" and either a ref/alt swap or a strand flip would be possible for a variant, the ref/alt swap has priority. This is only the case for A/T and C/G SNPs.

##### \-\-imputeInfo
Depending on how you want to further analyse the imputed data, you may want to look at the "\-\-imputeInfo" options. The default is "a" to output the allele dosages along with  imputed haplotypes. You can also add "g" for genotype dosages or "p" for genotype probabilities, or any combination of these three options if this suits your postprocessing tool better.
We recommend not to set all output types together because, first, all information can be easily obtained from the allele dosages and, second, the increased output has a significant impact on the size of the output files as well as an impact on the imputation speed.

##### \-\-imputeFilter
To reduce the size of the imputation output, you can filter all variants imputed below a certain R<sup>2</sup> threshold, which can be specified with this option.

##### -K
The "-K" option determines the number of best haplotypes selected from the reference for phasing. The default value is set to 10000, which is a good compromise between runtime and phasing quality in most cases. This value can be increased to achieve better phasing and thus better imputation quality, but with the disadvantage of significantly longer runtimes. If the argument is 0, K is chosen as maximum, which means that all haplotypes from the reference are used for phasing.
However, if K is greater than the number of haplotypes in the specified reference, all haplotypes will be used anyway and increasing K will have no effect.

##### \-\-maxChunkMemory
Since EagleImp performs automatic chunking of your data depending on memory requirements, you should provide the size of your system's RAM here. However, the default is set to 16 GB, and if the target data or the reference does not get too large, most imputations can be performed with much less memory needed.

##### \-\-skipPhasing
If the input data already contains phase information you can skip the phasing part with this switch to perform imputation only.

##### \-\-skipImputation
You can skip the imputation part with this switch to do a phasing of the input data only.

### Output files

##### \*.imputed.vcf.gz\/.bcf
The main output of the imputation contains both the phased target variants (overlap variants shared between target and reference) and the imputed variants from the reference. Each haplotype is provided with the corresponding allele code (0 or 1) and optionally with the dosage information according to the "\-\-imputeInfo" option, which is set to "a" (allele dosages) per default. (See the explanation above.)

The allele dosages of imputed variants are calculated from the imputation algorithm, while the dosages of the phased variants are calculated from the phase probability.

Furthermore, each variant contains information on the imputation quality score R<sup>2</sup>, allele frequency of the reference panel, allele frequency, minor allele frequency, allele count and allele number. Alleles and variant IDs are taken from the reference.

##### \*.phased.vcf.gz\/.bcf
Optional phasing output can be enabled with the "\-\-outputPhasedFile" switch. The file contains only variants from the phasing output. If this option is specified in combination with "\-\-outputUnphased", all variants from the target data set that were excluded from the analysis (for whatever reason) will also be included in the phasing output.

##### \*.varinfo
This file contains information for each target variant from the target data set, whether it was included or excluded from the analysis and for what reason.

##### \*.phased.confidences
This file contains the average phasing confidence for each sample.

##### Status files \/ Log file
The "\-\-stat" option allows you to select a file that will be continuously updated with the current progress status during analysis. If this option is enabled, an additional *\<statusfile\>.info* file is created to collect further information on the analysis (summarized in HTML-format from the normal program output). If necessary, warnings and errors from the running analysis are collected in a *\<statusfile\>.warning* and a *\<statusfile\>.error* file.

**NOTE:** All information provided in the status files is also presented in the normal program output. If you want to have a logfile, you can simply redirect the program output, e.g. directly or with the tool "tee", to split the output between the terminal and the log file (note that without the "2>&1" warning or error messages would not be logged!):
```
$ eagleimp <options> > mylogfile.log 2>&1
```
or
```
$ eagleimp <options> 2>&1 | tee mylogfile.log
```

### Notes

#### Phasing and imputation of chromosome X
The analysis of chromosome X is special when dealing with male (haploid) target samples. Usually, the pseudo-autosomal regions on chromosomes X and Y of male samples are presented together with the nonPAR region of chromosome X in a single chromosome X file. This results in diploid coding in the PAR regions but haploid coding in the nonPAR region. Consequently, phasing of male samples on chromosome X is only possible in the PAR regions, which is why **these regions should be handled independently of the nonPAR region.** For this reason, **EagleImp does not allow samples to be typed as diploid and haploid in the same file.** 
  
If you only have a single reference file ready for chromosome X that contains PAR and nonPAR regions together, you should split it with bcftools before creating a Qref file (the example includes region boundaries for GRCh37/hg19 build and assumes that the chromosome coding in your reference file is "X"):
```
$ bcftools view X.haplotypes.vcf.gz -r X:1-2699520 -Oz -o X_PAR1.haplotypes.vcf.gz
$ bcftools view X.haplotypes.vcf.gz -r X:2699521-154931043 -Oz -o X_nonPAR.haplotypes.vcf.gz
$ bcftools view X.haplotypes.vcf.gz -r X:154931044- -Oz -o X_PAR2.haplotypes.vcf.gz
```
Do the same for your target chromosome X file, or use the "launch_eagleimp" script, which will automatically split and then merge the target chromosome X file.
  
**NOTE for compatibility with the "launch-eagleimp" script:** Please name your Qref reference files according to the example, i.e. the chromosome identifier should be "X" or "23" suffixed with "_PAR1", "_nonPAR", "_PAR2", respectively. 
In any case, you do not need to create a separate genetic map for each region. It works well with one map for the complete chromosome X.

#### Imputation of chromosome Y
Since chromosome Y data is only available for male samples, it is expected that all data is haploid. Consequently, no phasing is required and EagleImp automatically skips the phasing step. As the genetic map is also only required in the phasing step, no genetic map of chromosome Y is required.
As with other chromosomes, you should note that your chromosome codings in the target and the reference must be the same when using a VCF reference. For a Qref reference, it does not matter whether you encode chromosome Y as "Y" or "24" with or without leading "chr".

#### How are missing haplotypes handled?
This depends on whether missings occur in the target data or in the reference. 
In the reference, missing haplotypes are randomly placed on the reference or alternative allele, distributed according to allele frequency. This is necessary because the PBWT data structure does not support missings. However, as long as the missing rate in the reference is low, this should have little effect.
For missings in the target, you can choose to impute them explicitly in the phasing process ("\-\-impMissing") or later during imputation. It makes more sense to do this in the imputation step (the default), since missings are **not** set to the reference allele and then used as anchors, as in the original PBWT tool. Instead, they are simply treated as if they were not present, and imputed like all other reference-only variants.

#### Phasing iterations
Currently, EagleImp performs only one phasing iteration by default. In most use cases, this is sufficient as long as the number of target samples does not significantly exceed the number of samples in the provided reference. (Our benchmarks have shown no quality improvement while the runtime increases linear with the number of iterations.) However, you can explicitly specify the number of iterations with the "-i" option.
In the first iteration, only the specified reference is used to determine the haplotype phases. From the second iteration on, all reference haplotypes are supplemented with the phased target haplotypes from the previous iteration, and the K best haplotypes are selected from the union of all reference and phased target haplotypes.
