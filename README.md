# Sarek Pipeline input

CLI for generating [Sarek pipeline](https://github.com/SciLifeLab/Sarek) input tsv file
for variant calling from a QBiC project stored in OpenBIS or ICGC data downloaded in AWS.

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

```bash
module load qbic/anaconda2/5.3.0
```

Execute without arguments or with -h to get an overview of the options:

```bash
$ python Sarek_pipeline_input.py -h

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

## Test

Under test/ you will find an example for Sarek input file preparation with empty fastq files.
The fastq files are in the directory test/fastq and the openbis databases in the directory test/openBIS/.
To produce the individual sarek input tsv files for each patient run:

```bash
$ python Sarek_pipeline_input.py -c Secondary_name -m -pT HCC -p test/fastas/ QMSHS test/openBIS/entity-browser-grid-sample.tsv test/openBIS/entity-browser-grid_experiments.tsv
Created merged data frame. Rows: 89 , Cols: 45
Eliminated empty columns. Rows: 89 , Cols: 36
Eliminated rows with no tissue annotation. Rows: 89 Cols: 36
Added boolean tumor annotation. Rows: 89 Cols: 38
Selected only DNA samples. Rows: 53 Cols: 38
Created an input file per patient.
```

Your FastQ files are now organized in folders following the structure:

```bash
/Entity(patient)/Biological_sample_code/(Tumor|Normal)/Test_sample_code/
```

The move command will warn you if there are extra DNA samples in your project's OpenBIS sample table that you didn't
download using qPostman. In this case 17 datasets where not downloaded.

The Sarek input files for each patient were created and stored in /test/fastqs. You can also have a look at your project
structure in the tree.txt file. You can have a look at the OpenBIS_merged.tsv file in case you need more information 
about other properties from the OpenBIS database for your samples (e.g. method used for exome capture, tissue where
normal samples come from).
