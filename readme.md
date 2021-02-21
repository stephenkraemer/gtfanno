# Introduction

This is yet another tool for annotating genomic regions. Why was that necessary?

1. This tool provides more detailed information than we found in existing tools.
2. We leverage all this information to resolve multiple potential annotations for the same regions into primary annotations using a set of precedence rules. This gives a single annotation per region, with a clear decision tree for each choice.
3. This tool provides a more flexibel definition of regions of interest with regard to TSS than we found in existing tools. For example, we like to partition the regions surrounding the TSS into a custom promoter region (say 5000 upstream and 1000 bp downstream of the TSS), and an adjacent more distant cis-regulatory region upstream of the promoter (maybe up to 45000 bp upstream of the promoter region). That's easy with gtfanno.

The tool provides the following basic annotations in this order of precedence (by default):
- Promoter
- 5'UTR, 3'UTR
- Intron, Exon
- Distant cis-regulatory region (see above)
- Intergenic

The basic annotation workflow is as follows
- For each region, find all possible annotations. An annotation is considered if the center of the region lies within the feature.
- For all annotations, various statistics are recorded, such as the number and percentages of overlapping bp or the distance to the TSS of the gene.
- For each region, define the primary feature class (e.g. Promoter, or Intron).
- Optionally, select the preferred annotation within the primary feature class according to various statistics. For example, a region may be annotated to the promoters of multiple genes. The top annotation could be chosen by the distance of the region center to the TSS, or by the percentage of overlapping bp.

This workflow produces the following output files:

Output files:
  - all_annotations (as BED3+ and pickled DataFrame)
    - contains all annotations
  - primary_annotations (BED3+ and pickled DataFrame):
    - contains only primary annotations for each region

Optional selection of the top annotation for each gene can be performed with simple dataframe operations based on the primary_annotations dataframe. See the usage example [here](./doc/usage.ipynb)

# Package maturity

This package is unreleased and unpublished software, but we use it often in in-house projects. We can not yet provide support for external users. Also, given the development status, we do change the API from time to time, without regard for health and safety of external users.


# Installation

There are no pypi or conda packages yet. The package can be installed directly from the repository:

```
pip install git+https://github.com/stephenkraemer/gtfanno.git
```

The package requires python >= 3.8 and pybedtools (and thus also: *bedtools*), see setup.py for all requirements. We recommend installing the package into a conda environment. Consider installing a tagged version for reproducibility. 

```
conda create -n gtfanno_env python=3.8 pybedtools
conda activate gtfanno_env
pip install git+https://github.com/stephenkraemer/gtfanno.git@0.2.0

```

# Supported operating systems

We only support Linux at the moment. The package may work on MacOS, but this is untested.

# Usage

The package provides a single public function:

```
ga.annotate(
    query_bed='/path/to/regions.bed',
    gtf_fp='/path/to/gencode.gtf',
    trunk_path='/path/to/output/files/basename',
    tmpdir='/path/to/tempdir',
    promoter=(-5000, 1000),
    distant_cis_regulatory_domain=(-50_000, -5000),
)
```

See the usage notebook (here)[./doc/usage.ipynb] for details.


