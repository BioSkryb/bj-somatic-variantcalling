
# BaseJumper BJ-SomaticVariantCalling

Pipeline processes WGS/Exome/Targeted sequencing data and performs comprehensive evaluation of single-cell libraries, and calls somatic SNP/Indel variants.

# Pipeline Overview
Pipeline requires normal bulk sample for each group of single-cell samples. If bulk sample is not provided, then the user must set "variant_workflow_type" to "pseudobulk" and pipeline will create a pseudobulk by subsampling uniformly across all single-cell samples. At minimum 3 single-cell samples are required per group to create normal psedobulk sample.

Following are the steps and tools that pipeline uses to perform the analyses:

- Map reads to reference genome using SENTIEON BWA MEM
- Remove duplicate reads using SENTIEON DRIVER LOCUSCOLLECTOR and SENTIEON DRIVER DEDUP
- Perform base quality score recalibration (BQSR) using SENTIEON DRIVER BQSR
- Perform variant calling with TNScope (defualt) or TNSEQ caller
- Perform variant annotation with SNPEFF, ClinVar, and dbSNP databases
- Evaluate metrics using SENTIEON DRIVER METRICS which includes Alignment, GC Bias, Insert Size, and Coverage metrics
- Aggregate the metrics across biosamples and tools to create overall pipeline statistics summary using MULTIQC

# Running Locally

Following are instructions for running BJ-SomaticVariantCalling in a local Ubuntu server

## Install Java

```
sudo apt-get install default-jdk

java -version
```

## Install AWS CLI

```
curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip"
unzip awscliv2.zip
sudo ./aws/install
```

## Install Nextflow

```
wget -qO- https://get.nextflow.io | bash
sudo mv nextflow /usr/local/bin/

```

## Install Docker

```
# Add Docker's official GPG key:
sudo apt-get update
sudo apt-get install ca-certificates curl
sudo install -m 0755 -d /etc/apt/keyrings
sudo curl -fsSL https://download.docker.com/linux/ubuntu/gpg -o /etc/apt/keyrings/docker.asc
sudo chmod a+r /etc/apt/keyrings/docker.asc

# Add the repository to Apt sources:
echo \
  "deb [arch=$(dpkg --print-architecture) signed-by=/etc/apt/keyrings/docker.asc] https://download.docker.com/linux/ubuntu \
  $(. /etc/os-release && echo "$VERSION_CODENAME") stable" | \
  sudo tee /etc/apt/sources.list.d/docker.list > /dev/null
sudo apt-get update

sudo apt-get install docker-ce docker-ce-cli containerd.io docker-buildx-plugin docker-compose-plugin

```

## Sentieon License Setup

The Sentieon license is a "localhost" license that starts a lightweight license server on the localhost. This type of license is very easy to use and get started with. However, because it can be used anywhere, we restrict this license to short-term testing/evaluation only. To use this type of license, you need to set the environment variable SENTIEON_LICENSE to point to the license file on the compute nodes:
```
export SENTIEON_LICENSE=</path/to/sentieon_eval.lic>
```
The license file should be saved at the base directory of the pipeline eg: `bj-somatic-variantcalling/sentieon_eval.lic`
All users will need to  [ submit helpdesk ticket](https://bioskryb.atlassian.net/servicedesk/customer/portal/3/group/14/create/156) to get an evaluation/full pass-through BioSkryb's Sentieon license.

## Resources Required

For running the pipeline, a typical dataset requires 8 CPU cores and 50 GB of memory. For larger datasets, you may need to increase the resources to 64 CPU cores and 120 GB of memory. You can specify these resources in the command as follows:
```
--max_cpus 8 --max_memory 50.GB
```
## Test Pipeline Execution

All pipeline resources are publically available at `s3://bioskryb-public-data/pipeline_resources` users need not have to download this, and will be downloaded during nextflow run.

**Command**

example-

** csv input **

```
git clone https://github.com/BioSkryb/bj-somatic-variantcalling.git
cd bj-somatic-variantcalling
nextflow run main.nf --input_csv $PWD/tests/data/inputs/input2.csv --publish_dir results/bj-somatic-variantcalling --max_cpus 8 --max_memory 50.GB
```
**Input Options**

The input for the pipeline can be passed via a input.csv with a meta data.

- **CSV Metadata Input**: The CSV file should have 6 columns: `biosampleName`, `read1`, `read2`, `groups`, `isbulk` and `bam`. 

Required Metadata:

- `biosampleName` column contains the name of the biosample. 
- `read1` and `read2` fields specify the paths to the input fastq files. Alternatively, the `bam` field provides the path to the bam file. Please note that either read1/read2 or bam must be provided. If the input is in bam format, read1 and read2 can be left blank. Conversely, if the input is in fastq format, the bam field can be left empty.  
- `groups` field contains the name of the sample group.  
- `isbulk` field is a boolean that indicates whether the sample is a bulk sample or a single/tumor sample.


For example:

Input.csv with fastq as input with 2 groups and both the groups have 1 bulk and 2 single/tumor sample each.

```
biosampleName,read1,read2,groups,isbulk,bam
chr22_testsample1_S1_L001,s3://bioskryb-public-data/pipeline_resources/dev-resources/local_test_files/chr22_testsample1_S1_L001_R1_001.fastq.gz,s3://bioskryb-public-data/pipeline_resources/dev-resources/local_test_files/chr22_testsample1_S1_L001_R2_001.fastq.gz,GROUP1,true,
chr22_testsample2_S1_L001,s3://bioskryb-public-data/pipeline_resources/dev-resources/local_test_files/chr22_testsample2_S1_L001_R1_001.fastq.gz,s3://bioskryb-public-data/pipeline_resources/dev-resources/local_test_files/chr22_testsample2_S1_L001_R2_001.fastq.gz,GROUP1,false,
chr22_testsample3_S1_L001,s3://bioskryb-public-data/pipeline_resources/dev-resources/local_test_files/chr22_testsample3_S1_L001_R1_001.fastq.gz,s3://bioskryb-public-data/pipeline_resources/dev-resources/local_test_files/chr22_testsample3_S1_L001_R2_001.fastq.gz,GROUP1,false,
chr22_testsample4_S1_L001,s3://bioskryb-public-data/pipeline_resources/dev-resources/local_test_files/chr22_testsample4_S1_L001_R1_001.fastq.gz,s3://bioskryb-public-data/pipeline_resources/dev-resources/local_test_files/chr22_testsample4_S1_L001_R2_001.fastq.gz,GROUP2,true,
chr22_testsample5_S1_L001,s3://bioskryb-public-data/pipeline_resources/dev-resources/local_test_files/chr22_testsample5_S1_L001_R1_001.fastq.gz,s3://bioskryb-public-data/pipeline_resources/dev-resources/local_test_files/chr22_testsample5_S1_L001_R2_001.fastq.gz,GROUP2,false,
chr22_testsample6_S1_L001,s3://bioskryb-public-data/pipeline_resources/dev-resources/local_test_files/chr22_testsample6_S1_L001_R1_001.fastq.gz,s3://bioskryb-public-data/pipeline_resources/dev-resources/local_test_files/chr22_testsample6_S1_L001_R2_001.fastq.gz,GROUP2,false,
```
Input.csv with fastq as input with 2 groups and both the groups have 3 single/tumor sample without a bulk sample. The parameter `--variant_workflow_type` should be set to `pseudobulk` for this run.

```
biosampleName,read1,read2,groups,isbulk,bam
chr22_testsample1_S1_L001,s3://bioskryb-public-data/pipeline_resources/dev-resources/local_test_files/chr22_testsample1_S1_L001_R1_001.fastq.gz,s3://bioskryb-public-data/pipeline_resources/dev-resources/local_test_files/chr22_testsample1_S1_L001_R2_001.fastq.gz,GROUP1,false,
chr22_testsample2_S1_L001,s3://bioskryb-public-data/pipeline_resources/dev-resources/local_test_files/chr22_testsample2_S1_L001_R1_001.fastq.gz,s3://bioskryb-public-data/pipeline_resources/dev-resources/local_test_files/chr22_testsample2_S1_L001_R2_001.fastq.gz,GROUP1,false,
chr22_testsample3_S1_L001,s3://bioskryb-public-data/pipeline_resources/dev-resources/local_test_files/chr22_testsample3_S1_L001_R1_001.fastq.gz,s3://bioskryb-public-data/pipeline_resources/dev-resources/local_test_files/chr22_testsample3_S1_L001_R2_001.fastq.gz,GROUP1,false,
chr22_testsample4_S1_L001,s3://bioskryb-public-data/pipeline_resources/dev-resources/local_test_files/chr22_testsample4_S1_L001_R1_001.fastq.gz,s3://bioskryb-public-data/pipeline_resources/dev-resources/local_test_files/chr22_testsample4_S1_L001_R2_001.fastq.gz,GROUP2,false,
chr22_testsample5_S1_L001,s3://bioskryb-public-data/pipeline_resources/dev-resources/local_test_files/chr22_testsample5_S1_L001_R1_001.fastq.gz,s3://bioskryb-public-data/pipeline_resources/dev-resources/local_test_files/chr22_testsample5_S1_L001_R2_001.fastq.gz,GROUP2,false,
chr22_testsample6_S1_L001,s3://bioskryb-public-data/pipeline_resources/dev-resources/local_test_files/chr22_testsample6_S1_L001_R1_001.fastq.gz,s3://bioskryb-public-data/pipeline_resources/dev-resources/local_test_files/chr22_testsample6_S1_L001_R2_001.fastq.gz,GROUP2,false,
```
Input.csv with bam as input with 2 groups and both the groups have 1 bulk and 2 single/tumor sample each. Please make sure to set the parameter `is_bam` to true when passing bam files as input.

```
biosampleName,read1,read2,groups,isbulk,bam
chr22_testsample1_S1_L001,,,GROUP1,true,s3://bioskryb-public-data/pipeline_resources/dev-resources/local_test_files/chr22_testsample1_S1_L001.bam
chr22_testsample2_S1_L001,,,GROUP1,false,s3://bioskryb-public-data/pipeline_resources/dev-resources/local_test_files/chr22_testsample2_S1_L001.bam
chr22_testsample3_S1_L001,,,GROUP1,false,s3://bioskryb-public-data/pipeline_resources/dev-resources/local_test_files/chr22_testsample3_S1_L001.bam
chr22_testsample4_S1_L001,,,GROUP2,true,s3://bioskryb-public-data/pipeline_resources/dev-resources/local_test_files/chr22_testsample4_S1_L001.bam
chr22_testsample5_S1_L001,,,GROUP2,false,s3://bioskryb-public-data/pipeline_resources/dev-resources/local_test_files/chr22_testsample5_S1_L001.bam
chr22_testsample6_S1_L001,,,GROUP2,false,s3://bioskryb-public-data/pipeline_resources/dev-resources/local_test_files/chr22_testsample6_S1_L001.bam
```
Input.csv with bam as input with 2 groups and both the groups have 3 single/tumor sample without a bulk sample. The parameters `is_bam` and  `--variant_workflow_type` should be set to `pseudobulk` for this run.

```
biosampleName,read1,read2,groups,isbulk,bam
chr22_testsample1_S1_L001,,,GROUP1,false,s3://bioskryb-public-data/pipeline_resources/dev-resources/local_test_files/chr22_testsample1_S1_L001.bam
chr22_testsample2_S1_L001,,,GROUP1,false,s3://bioskryb-public-data/pipeline_resources/dev-resources/local_test_files/chr22_testsample2_S1_L001.bam
chr22_testsample3_S1_L001,,,GROUP1,false,s3://bioskryb-public-data/pipeline_resources/dev-resources/local_test_files/chr22_testsample3_S1_L001.bam
chr22_testsample4_S1_L001,,,GROUP2,false,s3://bioskryb-public-data/pipeline_resources/dev-resources/local_test_files/chr22_testsample4_S1_L001.bam
chr22_testsample5_S1_L001,,,GROUP2,false,s3://bioskryb-public-data/pipeline_resources/dev-resources/local_test_files/chr22_testsample5_S1_L001.bam
chr22_testsample6_S1_L001,,,GROUP2,false,s3://bioskryb-public-data/pipeline_resources/dev-resources/local_test_files/chr22_testsample6_S1_L001.bam
```


**Optional Parameters**: 

- **variant_workflow_type**: This parameter sets the appropriate workflow to run based on the provided samples:
                - match_normal: if matched normal is provided for the set of single/tumor sample.
                - panel_normal: if there are panel of normal samples provided along with set of single cell/tumor samples.
                - pseudobulk: if no normal samples are provided.The pseudo bulk workflow is an optional feature that generates a pseudo bulk fastq file when all provided inputs are single/tumor samples
                - somatic_heurestic_filter: Uses dnascope to make variant calls and then uses a heuristic filter to keep somatic variants and filter out germline. It doesn't require bulk samples and creates a phylogenetic tree.

- **somatic_variant_caller**: There's another optional parameter available, `--somatic_variant_caller`, which allows users to select the Variant Caller. The options are `tnscope` or `tnseq`, with tnscope being the default choice.

- **panel_of_normal_vcf**: This parameter allows you to specify a panel of VCF files when running the `panel_normal` workflow. If this parameter is not provided, a panel of normal VCF will be created by default.

**Optional Modules**


This pipeline includes optional modules. You can choose to include or exclude these modules by adjusting the following parameters:

- `--skip_variant_annotation`: Set this to `true` to exclude the Variant Annotation module. By default, it is set to `true`.
- `--skip_sigprofile`: Set this to `false` to enable the Mutational Signature module. By default, it is set to `false`.


**Outputs**


The pipeline saves its output files in the designated "publish_dir" directory. 

- secondary_analysis/ : alignment/ The bam files after alignment are stored in the "secondary_analyses/alignment/" subdirectory  
                        metrics/   The metrics files are saved in the "secondary_analyses/metrics/<sample_name>_<group>/" subdirectory.  

- variant_calls_<tnscope / tnseq> : The vcf files after the variant calling are saved in the "variant_calls_<tnscope/tnseq>/<sample_name>_<group>/" subdirectory.  

- multiqc/ : This section includes output files containing metrics from various tools to create a MultiQC report. 

**command options**

`nextflow run main.nf --help`

```
BJ-SomaticVariantCalling

Usage:
    nextflow run main.nf [options]

Script Options: see nextflow.config

    [required]

    --input_csv                 FILE    Path to input csv file


    --publish_dir               DIR     Path of the output directory

    --genome                    STR     Reference genome to use. Available options - GRCh37, GRCh38
                                        DEFAULT: GRCh38

    [optional]

    --variant_workflow_type     STR    This parameter sets the appropriate workflow to run based on the provided samples:
                                            - match_normal: if one normal sample is provided for the set of single/tumor samples.
                                            - panel_normal: if more than one normal sample is provided for the set of single/tumor samples.
                                            - pseudobulk: if no normal samples are provided.
                                            - somatic_heurestic_filter: Uses dnascope to make variant calls and then uses a heuristic filter to keep somatic variants and filter out germline. It doesn't require bulk samples.

    --panel_of_normal_vcf      FILE    Path to panel of normal vcf file

    --somatic_variant_caller    STR     To select the Variant Caller. The options are tnscope or tnseq.
                                        DEFAULT: tnscope

    --multiqc_config            DIR     Path to multiqc config
                                        DEFAULT: bj-somatic-variantcalling/assets/multiqc

    --skip_variant_annotation   BOOL    Whether to skip variant annotation
                                        DEFAULT: null
                                        
    --skip_sigprofile           BOOL    Whether to skip Mutational Signature
                                        DEFAULT: false

    --help                      BOOL    Display help message

```
**Tool versions**

- `Sentieon: 202308.01`
- `VEP: 111.0`
- `GATK: 4.1.3.0`
- `BCFtools: 1.14`
- `VariantQC: 1.20`

**nf-test**

The BioSkryb bj-somatic-variantcalling nextflow pipeline run is tested using the nf-test framework.

Installation:

nf-test has the same requirements as Nextflow and can be used on POSIX compatible systems like Linux or OS X. You can install nf-test using the following command:
```
wget -qO- https://code.askimed.com/install/nf-test | bash
sudo mv nf-test /usr/local/bin/
```
It will create the nf-test executable file in the current directory. Optionally, move the nf-test file to a directory accessible by your $PATH variable.

Usage:

```
nf-test test
```

The nf-test for this repository is saved at tests/ folder.

```
    test("somatic_variant_calling_test") {

        when {
            params {
                publish_dir = "${outputDir}/results"
                input_csv   = "$baseDir/tests/data/inputs/input2.csv"
                timestamp = "test"
                variant_workflow_type = "pseudobulk"
                is_bam = "true"
                architecture = "x86"
            }
        }

        then {
            assertAll(
                // Check if the workflow was successful
                { assert workflow.success },

                // Verify existence of the multiqc report HTML file
                {assert new File("${outputDir}/results_test/multiqc/multiqc_report.html").exists()},

                // Check for a match in the all metrics MQC text file
                {assert snapshot (path("${outputDir}/results_test/secondary_analyses/metrics/nf-wgs-pipeline_all_metrics_mqc.txt")).match("all_metrics_mqc")},

                // Check for a match in the selected metrics MQC text file
                {assert snapshot (path("${outputDir}/results_test/secondary_analyses/metrics/nf-wgs-pipeline_selected_metrics_mqc.txt")).match("selected_metrics_mqc")},

                // Verify existence of the bam file
                {assert new File("${outputDir}/results_test/secondary_analyses/alignment/GROUP1.bam.bai").exists()},

                // Check for a match in the pseudobulk csv file
                {assert snapshot (path("${outputDir}/results_test/pseudobulk/pseudobulk_input.csv")).match("pseudo_bulk_input")},

                // Verify existence of the tnscope out file
                {assert new File("${outputDir}/results_test/variant_calls_tnscope/CD-Cas2-S12_S1_L001_GROUP1/CD-Cas2-S12_S1_L001_tnscope.vcf.gz.tbi").exists()}

            )
        }

    }
```

# Need Help?
If you need any help, please [submit a helpdesk ticket](https://bioskryb.atlassian.net/servicedesk/customer/portal/3/group/14/create/156).

# References
For more information, you can refer to the following publications:

- Chung, C., Yang, X., Hevner, R. F., Kennedy, K., Vong, K. I., Liu, Y., Patel, A., Nedunuri, R., 
Barton, S. T., Noel, G., Barrows, C., Stanley, V., Mittal, S., Breuss, M. W., Schlachetzki,
J. C. M., Kingsmore, S. F., & Gleeson, J. G. (2024). Cell-type-resolved mosaicism 
reveals clonal dynamics of the human forebrain. Nature, 629(8011), 384–392. [https://doi.org/10.1038/s41586-024-07292-5](https://doi.org/10.1038/s41586-024-07292-5)

- Zhao, Y., Luquette, L. J., Veit, A. D., Wang, X., Xi, R., Viswanadham, V. V, Shao, D. D., 
Walsh, C. A., Yang, H. W., Johnson, M. D., & Park, P. J. (2024).
High-resolution 
detection of copy number alterations in single cells with HiScanner. BioRxiv, 
2024.04.26.587806. [https://www.biorxiv.org/content/10.1101/2024.04.26.587806v1.full](https://www.biorxiv.org/content/10.1101/2024.04.26.587806v1.full)

- Zawistowski, J. S., Salas-González, I., Morozova, T. V, Blackinton, J. G., Tate, T., Arvapalli, 
D., Velivela, S., Harton, G. L., Marks, J. R., Hwang, E. S., Weigman, V. J., & West, J. A. 
A. (n.d.). Unifying genomics and transcriptomics in single cells with ResolveOME amplification
chemistry to illuminate oncogenic and drug resistance mechanisms. [https://www.biorxiv.org/content/10.1101/2022.04.29.489440v1.full](https://www.biorxiv.org/content/10.1101/2022.04.29.489440v1.full)

NOTE: Several studies have utilized BaseJumper pipelines as part of the standard 
quality control processes implemented through ResolveServices<sup>SM</sup>. While these 
pipelines may not be explicitly cited, they are integral to the methodologies 
described.
