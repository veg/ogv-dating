# Timing of QVOA strains using RNA pre-treatment data

<a href="https://zenodo.org/badge/latestdoi/181484600"><img src="https://zenodo.org/badge/181484600.svg" alt="DOI"></a>
> Updated Nov 2022 to provide compatibility with HyPhy 2.5.x and support additional options, like not compressing identical sequences
> Updated Mar 2024 to provide compatibility with HyPhy 2.5.5x, better handle alignment when none of the sequences are in frame (use longest ORF), and branch length re-estimation using HyPhy.

## Overview

This pipeline **times** outgrowth virus (OGV) strains from a single host using longtitudinally sampled RNA data to provide temporal information. Four different approaches are used to time strains with unobserved dates (OGV or DNA strains). The phylogenetic tree for each genomic segment is rooted to maximize the root-to-tip to sampling time correlation coefficient first.

### Patristic distance
Given a tree with OGV and RNA strains, this approach simply finds the nearest RNA neighbor (in [patristic](https://link.springer.com/referenceworkentry/10.1007/978-1-4020-6754-9_12406), i.e., path length sense) for each OGV sequence. To assess reliability of this assignment, we also identify all **K≥1** sequences within a factor of `2x` of this smallest distance, and extract the consensus vote to assign the date. The **support** value if the proportion of **K** neighboring sequences that come from the same timepoint. In the example below, there are **9** RNA neighbors to the blue QVOA sequence. The closest one is `CAP288_4260_184WPI_ENV_C1C2_NGS_020_0_006` and the other two are with 2x distance. The blue sequence will be assigned to the `184WPI` timepoint (WPI = weeks post infection), because most of the neighbors are from that timepoint, and the support value will be `8/9`, because one of teh neighbors has a timestamp of `174WPI`

![img1](img/1.png)

### Clade support
Here, we start with a QVOA sequence, and go up the tree until we encounter the first internal node that has **both** QVOA and RNA descendants and that has a minimum of `0.90` bootstrap support. We then assign a timestamp and support to the QVOA sequence using the same majority rule as for the patristic distance for all the nodes in the clade. In the example below, to date `CAP288_NEF_1_W11_QVOA`, we walk up the tree to the first supported node (blue) and look at its RNA descendants, of which there are **6** AND 5/6 come from `184WPI`, yielding the assignment of `184WPI` for the QVOA sequence and support of 5/6.

![img1](img/2.png)

### Placement

Phylogenetic placement using a modified [IgSCUEAL](https://royalsocietypublishing.org/doi/full/10.1098/rstb.2014.0240) tool. The support value is the model weight (in [0-1]) assigned to the branches with teh corresponding label

### Regression

The linear regression molecular clock method, essentially the same as in [the Jones et al](https://www.pnas.org/content/115/38/E8958) paper.

<p>

## Implementation 

### Dependancies

1. `HyPhy ≥v2.5.55` -- or the `develop` branch from [github.com/veg/hyphy](github.com/veg/hyphy)
2. `MAFFT` for MSA generation
3. `FastTree` for phylogeny inference
4. `phylotree.js` (an `npm` package) for tree processing and visualization [note; there is a custom JS script for performing dating steps that is included in `scripts`]
5. `snakemake` for workflow management
6. `node` for executing JavaScript
7. `BioPython` 

See notes on configuration below and please run 
```
npm install 
```
from the base directory 


#### Workflow

The workflow for sample classification is implemented as a [Snakemake](https://snakemake.readthedocs.io/en/stable/) workflow, encoded in the `Snakemake` file, and can be graphically summarized as follows (for a specific sample)...

![img1](img/3.png)

The pipeline requires an input `samples.json` file, with pathnames of files to process specified relative to the `data` directory, e.g. 

```
{
    'samples': [
        "CAP188/CAP188_ENV_2_all_hap.fasta",
        "CAP188/CAP188_ENV_3_all_hap.fasta",
        "CAP188/CAP188_ENV_4_all_hap.fasta",
        "CAP188/CAP188_NEF_1_all_hap.fasta",
        "CAP206/CAP206_ENV_2_all_hap.fasta",
        "CAP206/CAP206_ENV_4_all_hap.fasta",
        "CAP206/CAP206_GAG_1_all_hap.fasta",
        "CAP217/CAP217_ENV_2_all_hap.fasta",
        "CAP217/CAP217_ENV_3_all_hap.fasta",
        "CAP217/CAP217_ENV_4_all_hap.fasta",
        "CAP257/CAP257_ENV_2_all_hap.fasta",
        "CAP257/CAP257_ENV_3_all_hap.fasta",
        "CAP257/CAP257_ENV_4_all_hap.fasta",
```

The pipeline is executed by calling `snakemake`, e.g. 

```
$snakemake --cores 4 
```

After the pipeline has completed, you need to run 

```
$python3 scripts/result-summary.py -d results/dating/ -j data/conversion.json -o results/summary.csv
```

to generate the overall summary CSV file.

---

The high level directory structure contains the following components

```
├── HBL							[HyPhy scripts]
│   ├── data_filter-2.bf
│   ├── data_filter.bf
│   ├── lib
│   └── scripts
├── LICENSE
├── package.json					[npm package manifest]
├── README.md						[This file]
├── Snakefile						[Snakemake configuration file]
├── data							[unprocessed FASTA files]
│   ├── CAP188
│   ├── CAP206
│   ├── CAP217
│   ├── CAP257
│   ├── CAP268
│   ├── CAP280
│   ├── CAP287
│   ├── CAP288
│   ├── CAP302
│   ├── CAP316
│   ├── CAP336
│   ├── CAP372
│   └── conversion.json	  [Color codes and WPI to pre-ART conversions]
├── img					  [Images for this file]
│   ├── 1.png
│   ├── 2.png
│   └── 3.png
├── results              [Results generated by the pipeline]
│   ├── alignments	     [Alignments]
│   ├── dating		     [QVOA/OGV dating results]
│   └── trees            [Trees inferred by the pipeline]
├── samples.json         [Manifest of files to be processed by the pipeline]
└── scripts
    ├── compute-distance.js 	[Classifier/timing script]
    └── result-summary.py		[Penultimate summary script]
```    


Within specific resutls directories, the following files are created **for each input FASTA**

#### `alignments` directory
```
CAP316/
├── CAP316_ENV_3_all_hap.fasta_combined.msa           
    [RNA+QVOA nucleotide in-frame sequences, aligned]
├── CAP316_ENV_3_all_hap.fasta_combined_nuc.fas
    [RNA+QVOA nucleotide in-frame sequences, not aligned]
├── CAP316_ENV_3_all_hap.fasta_combined_protein.fas
    [Translated RNA+QVOA nucleotide in-frame sequences, not aligned]
├── CAP316_ENV_3_all_hap.fasta_combined_protein.msa
    [Translated RNA+QVOA nucleotide in-frame sequences, aligned]
├── CAP316_ENV_3_all_hap.fasta_nuc.fas
    [In frame RNA nucleotide in-frame sequences, not aligned]
├── CAP316_ENV_3_all_hap.fasta_protein.fas
    [Translated RNA nucleotide in-frame sequences, not aligned]
├── CAP316_ENV_3_all_hap.fasta_qvoa.fas
    [QVOA/OGV in-frame sequences, not aligned]
├── CAP316_ENV_3_all_hap.fasta_rna.msa
    [RNA+QVOA nucleotide in-frame sequences not aligned]
...

```

#### `trees` directory

```
CAP316/
├── CAP316_ENV_3_all_hap.fasta_combined_unscaled.nwk
    [Pseudo-ML FastTree topology for QVOA + RNA sequences]
├── CAP316_ENV_3_all_hap.fasta_rna_unscaled.nwk
    [Pseudo-ML FastTree topology for RNA sequences]
├── CAP316_ENV_3_all_hap.fasta_combined.nwk
    [FastTree topology with branch lengths retuned by HyPhy]
├── CAP316_ENV_3_all_hap.fasta_rna.nwk
    [FastTree topology with branch lengths retuned by HyPhy]
```

#### `dating` directory

```
CAP336/
├── CAP336_ENV_4_all_hap.fasta.scueal
	[alignment+tree file used by HyPhy for phylogenetic placement]
├── CAP336_ENV_4_all_hap.fasta.scueal.hyphy
	[auxiliary file generated by HyPhy during phylogenetic placement]
├── CAP336_ENV_4_all_hap.fasta.scueal.hyphy-data
	[auxiliary file generated by HyPhy during phylogenetic placement]
├── CAP336_ENV_4_all_hap.fasta_classification.csv
	[Summary classification (see below)]
├── CAP336_ENV_4_all_hap.fasta_full_placement.tsv
	[phylogenetic placement extended results]
├── CAP336_ENV_4_all_hap.fasta_placement.tsv
	[phylogenetic placement brief results]
```

#### Alignment preparation

For each individual and each genomic region, we performed the following steps to generate in-frame codon alignments for subsequent analysis. Sequence names in the input files are expected to encode their origin (RNA vs QVOA): any sequence whose name matched the regexp `(QVOA)|(OGV)|(DNA)|(xxxx)|(BAR[0-9]+)` are interpreted as outgrowth virus (see https://github.com/veg/ogv-dating/blob/128af4b21e20d28e76c5848eeac30e06a6be6d98/HBL/data_filter.bf#L41, for example).

```
CAP188_xxxx_000WPI_ENV_C2C3_OGV_B-W36
                            --- [OGV]
CAP188_4260_185WPI_ENV_C2C3_NGS_046_0.002
            ------ [Time stamp in weeks post-infection]
```
1. Translate each of the **RNA** nucleotide sequences into protein sequences using 3 forward reading frames. Keep track those that are in frame. _Note : for some of the alignments **none** of the sequences were in frame, typically because the reads were **stitched** from discontinuous genomic regions [e.g. due to variable loop removal]. For such alignments, the code will find the longest ORF to use as a reference, meaning that only partial sequence data may be retained_


2. Reads that were **not** in frame, and QVOA reads were codon aligned to **each** of the in-frame reads using a codon-aware Smith-Waterman algorithm implemented in `HyPhy`, picking the best-scoring pair to correct frame-shifts.

3. The resulting in-frame sequences were translated to amino-acids, mapping all ambiguous codons to `?` (e.g. `A-A`).

4. A multiple sequence alignment (MSA) was generated from translated protein sequences using `mafft`

5. Original codon sequences are mapped to aligned protein sequences using a script in `HyPhy`. Identical sequences (which do not contribute any information to phylogenetic analyses) are removed at this stage.

6. ML phylognetic trees are reconstructed using `FastTree` (option) for RNA sequences **only** (these are later used for phylogenetic placement) and for RNA+QVOA. These topologies are then fed into `hyphy` to re-estimate the branch lengths.

7. The pipeline will **compress** identical RNA sequences, i.e. only one copy will be retained. The copy number will (including single copies, i.e. `1`) will be appended to sequence names. It will NOT compress QVOA sequences, defined as those with `QVOA` or `OGV` in their name.

> **Important:** adjust the paths to various executables, like `hyphy` and `FastTree`, in the `Snakefile` before running the pipeline

```
callables = {
    'hyphy' : '/usr/local/bin/hyphy',
    'mafft' : '/usr/local/bin/mafft --auto ',
    'raxml' : '/usr/local/bin/raxml-ng',
    'fasttree' : '/usr/local/bin/FastTree',
    'classifier' : 'scripts/compute-distance.js'
}
```