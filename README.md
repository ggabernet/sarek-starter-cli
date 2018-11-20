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
module load qbic/anaconda2/5.3.0
```

Execute without arguments or with -h to get an overview of the options:
```
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

```
$ python Sarek_pipeline_input.py -c Secondary_name -m -pT HCC -p test/fastas/ QMSHS test/openBIS/entity-browser-grid-sample.tsv test/openBIS/entity-browser-grid_experiments.tsv
Created merged data frame. Rows: 89 , Cols: 45
Eliminated empty columns. Rows: 89 , Cols: 36
Eliminated rows with no tissue annotation. Rows: 89 Cols: 36
Added boolean tumor annotation. Rows: 89 Cols: 38
Selected only DNA samples. Rows: 53 Cols: 38
mv: rename I17D005b01_01* to QMSHSENTITY-11/QMSHS869AI/Normal/QMSHS271CB/I17D005b01_01*: No such file or directory
mv: rename I17D005b02_01* to QMSHSENTITY-11/QMSHS870AL/Normal/QMSHS272CJ/I17D005b02_01*: No such file or directory
mv: rename I17D005b03_01* to QMSHSENTITY-11/QMSHS896AF/Normal/QMSHS296CV/I17D005b03_01*: No such file or directory
mv: rename I17D005b05_01* to QMSHSENTITY-12/QMSHS950AD/Normal/QMSHS273CR/I17D005b05_01*: No such file or directory
mv: rename I17D005b06_01* to QMSHSENTITY-12/QMSHS951AL/Normal/QMSHS297C5/I17D005b06_01*: No such file or directory
mv: rename I17D005b07_01* to QMSHSENTITY-13/QMSHS081BV/Normal/QMSHS298CD/I17D005b07_01*: No such file or directory
mv: rename I17D005b08_01* to QMSHSENTITY-13/QMSHS082B5/Normal/QMSHS274C1/I17D005b08_01*: No such file or directory
mv: rename I17D005b11_01* to QMSHSENTITY-16/QMSHS324BV/Normal/QMSHS300CS/I17D005b11_01*: No such file or directory
mv: rename I17R018Db02_01* to QMSHSENTITY-17/QMSHS517BQ/Normal/QMSHS285CG/I17R018Db02_01*: No such file or directory
mv: rename I17R018Da02_02* to QMSHSENTITY-17/QMSHS558BS/Tumor/QMSHS279C7/I17R018Da02_02*: No such file or directory
mv: rename I17R018Da03_01* to QMSHSENTITY-17/QMSHS583B9/Tumor/QMSHS484CK/I17R018Da03_01*: No such file or directory
mv: rename I17R018Dd02_01* to QMSHSENTITY-18/QMSHS599BU/Normal/QMSHS481CU/I17R018Dd02_01*: No such file or directory
mv: rename I17R018Dd01_01* to QMSHSENTITY-18/QMSHS606BH/Tumor/QMSHS478CB/I17R018Dd01_01*: No such file or directory
mv: rename I17R018Df03_01* to QMSHSENTITY-22/QMSHS868BJ/Normal/QMSHS487CA/I17R018Df03_01*: No such file or directory
mv: rename I16D016a01_01* to QMSHSENTITY-5/QMSHS370AP/Tumor/QMSHS270C3/I16D016a01_01*: No such file or directory
mv: rename I16D016a02_01* to QMSHSENTITY-5/QMSHS379AT/Normal/QMSHS063BX/I16D016a02_01*: No such file or directory
mv: rename GS160323_02* to QMSHSENTITY-7/QMSHS497AX/Normal/QMSHS034BK/GS160323_02*: No such file or directory
Created an input file per patient.

```
Your FastQ files are now organized in folders following the structure:

```
/Entity(patient)/Biological_sample_code/(Tumor|Normal)/Test_sample_code/
```

The move command will warn you if there are extra DNA samples in your project's OpenBIS sample table that you didn't
download using qPostman. In this case 17 datasets where not downloaded.

The Sarek input files for each patient were created and stored in /test/fastqs. You can also have a look at your project
structure in the tree.txt file. You can have a look at the OpenBIS_merged.tsv file in case you need more information 
about other properties from the OpenBIS database for your samples (e.g. method used for exome capture, tissue where
normal samples come from).