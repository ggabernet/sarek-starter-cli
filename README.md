# Sarek Pipeline input

CLI for generating [Sarek pipeline](https://github.com/SciLifeLab/Sarek) input tsv file
for variant calling from a QBiC project stored in OpenBIS.

## Pre-requisites

* Download the fastq files for all the samples that you want to call variants for, 
for example using [qPostman](https://github.com/qbicsoftware/postman-cli).

* Download the all samples and all experiments tsv table for your project in OpenBIS.
  * Samples tsv: in OpenBIS Browse --> samples. Needed columns: Code, Parents, 
  Sample Type, Sample type, Project, Experiment, Additional information, 
  Secondary name, Primary tissue/body fluid. Export as tsv.
  * Experiments tsv: in OpenBIS Browse --> experiments. Look for your project. Needed
  columns: Experiment, Project, Code, Additional information.
  
## Usage
Load the conda package or your favourite python package manager (python 2.7 and pandas package is required). E.g. in new CFC cluster:

```
module load qbic/anaconda2/2.1.0
```

Execute without arguments or with -h to get an overview of the options:
```
~$ python Sarek_pipeline_input.py -h

usage: Sarek_pipeline_input.py [-h] [-p PATH] [-c {Test,Secondary_name}]
                               [-pR1 PATTERN_R1] [-pR2 PATTERN_R2]
                               [-pL PATTERN_LANE] [-m]
                               project sample_tsv experiment_tsv

positional arguments:
  project               QBiC project code for which VC should be calculated.
  sample_tsv            Path to the sample table tsv file extracted from
                        OpenBIS.
  experiment_tsv        Path to the experiment table tsv file extracted from
                        OpenBIS.

optional arguments:
  -h, --help            show this help message and exit
  -p PATH, --path PATH  Path to folder with fastq files.
  -c {Test,Secondary_name}, --contains {Test,Secondary_name}
                        String of the identifier that is contained in the
                        fastq filename. 'Test' stands for QBiC test sample
                        code. 'Secondary_name' stands for NGS sample secondary
                        name (default)(usu. Genetics ID).
  -pR1 PATTERN_R1, --pattern_R1 PATTERN_R1
                        Regex to look for at fastq filename and identify 1st
                        fastq of a pair.
  -pR2 PATTERN_R2, --pattern_R2 PATTERN_R2
                        Regex to look for at fastq filename and identify 2nd
                        fastq of a pair.
  -pL PATTERN_LANE, --pattern_lane PATTERN_LANE
                        Regex to look for at fastqfilename to
                        identifysequencing lane.
  -m, --multiple        Create a separate input file for each entity/patient.

```
