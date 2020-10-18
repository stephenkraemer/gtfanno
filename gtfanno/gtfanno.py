"""Annotate regions with gene annotations from a GTF file

Currently annotated features
- promoter
- exon, intron
- 5'-UTR, 3'UTR
- intergenic
- distant cis-regulatory domain (DCRD)

Input file restrictions
- query_bed: (this will be improved in the future)
  - no chromosome prefix allowed
  - must be sorted -k1,1 -k2,2n -k3,3n


Output files:
  - all_annotations (BED3+annocols, pickled DataFrame)
    - Saved to output_trunk + '_all-annotations.{bed,p}'
    - contains all annotations, classified as 'primary' or 'secondary' (column feature_rank).
  - primary annotations (BED3+annocols, pickled DataFrame):
    - saved to output_trunk + '_primary-annotations.{bed,p}'
    - contains only primary annotations for each region. Each region can
    only have one primary annotation per gene. A region may have
    multiple primary annotations if the same feature is annotated for
    several genes.

Annotation workflow:
  1. Assign each region all possible annotations
     - a query region is assigned to a promoter, exon, intron, UTR or DCRD if the center
       of the region lies within the feature
  2. Classify the annotations into primary and secondary
     - for each region, the annotated feature with the highest precedence is
       chosen as primary annotated feature. The underlying precedence list is:
       - Promoter
       - 5'-UTR
       - 3'-UTR
       - exon
       - intron
       - DCRD
     - the result of this comprehensive annotation step is saved in the all_annotations output
  2.

Notes on this workflow
- antisense dmrs
- alternative transcripts

Example call:

    annotate(query_bed='/icgc/dkfzlsdf/analysis/hs_ontogeny/results/wgbs/cohort_results/dmrs/pop-group-dmrs/full-hierarchy-dmrs_v1.bed',
         gtf_fp=("/icgc/dkfzlsdf/analysis/hs_ontogeny/databases/gene_annotations"
                 "/gencode.vM19.annotation.no-prefix.gtf"),
         trunk_path='/home/kraemers/temp/test',
         tmpdir='/icgc/dkfzlsdf/analysis/B080/kraemers/temp',
         promoter = (-5000, 1000),
         distant_cis_regulatory_domain = (-50000, -5000),
         )

Papercuts
- column names: feat_class, feature_rank, ...

Possible improvements
- currently, annotation columns in the query BED are discarded. Add an
  option to include them in the output file.
- allow the user to choose the transcript parts to be annotated. Currently hardcoded
  to exon, intron, 5'-UTR, 3'-UTR (not start_codon, ...)
- the feature precedence is hardcoded. Turn into parameter.
- when several annotations for a feature with the same precedence exist, the distance
  from the query center to the feature is used to select the best candidate. Allow using
  other statistics, e.g. percent_overlap, and make available as parameter. Different
  other statistics are already computed, so it is just a matter of using them.
- add an option to prefer the principal transcript when the same feature
  is annotated for several transcripts of a gene

Assumptions about the gene annotation format:
  - see the docstring of annotate_transcript_parts for assumptions and
    workflow choices for the annotation of transcript parts (exons, introns, ...)
"""

# %%
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Tuple

import numpy as np
import pandas as pd
from pandas.api.types import CategoricalDtype
from pybedtools import BedTool

# %%


def annotate(
    query_bed,
    gtf_fp,
    trunk_path,
    tmpdir,
    promoter=(-5000, 1000),
    distant_cis_regulatory_domain=(-50000, -5000),
):
    """Annotate query regions with feature from a GTF or GFF file

    Args:
        query_bed: Path to BED file with query regions. Expects only standard BED3 conventions,
        ie first three columns according to BED format. Header line is allowed if commented with '#'.
        Additional columns are ignored
        gtf_fp: Path to gtf or gff file
        trunk_path: output file path. Besides BED format, the output table
            will also be saved as pickled pandas.DataFrame
        tmpdir: for creation of temporary files (will be cleaned up
            automatically)
        promoter: (upper boundary, lower boundary) - both given as
            distance from the TSS. Negative numbers indicate upstream,
            positive number downstream boundaries.
        distant_cis_regulatory_domain: Same as promoter.
            Typically, the distant_cis_regulatory_domain is chosen to
            be directly adjacent to the proximal promoter.

    Results:
        Annotated regions in BED and pickled pandas.DataFrame format
    """

    assert Path(
        trunk_path
    ).parent.exists(), "The output directory does not exist, abort..."

    # %% Setup
    # --------------------------------------------------------------------------

    # Typical error is to specify upstream bp with negative number
    assert distant_cis_regulatory_domain[0] < 0, promoter[0] < 0
    # The upstream area of the DCRD should be larger than the one of the promoter
    assert distant_cis_regulatory_domain[0] < promoter[0]

    output_cols = [
        "Chromosome",
        "Start",
        "End",
        "gtfanno_uid",
        "center",
        "feat_class",
        "perc_feature",
        "perc_region",
        "distance",
        "has_center",
        "gene_name",
        "gene_id",
        "transcript_id",
        "appris_principal_score",
        "feat_chrom",
        "feat_start",
        "feat_end",
        "feat_center",
        "feat_strand",
    ]

    tmpdir_obj = TemporaryDirectory(dir=tmpdir)
    tmpdir_path = Path(tmpdir_obj.name)

    print("Loading data")
    # query df columns: Chromosome, Start, End, gtfanno_uid
    # the gtfanno_uid is kept across all operations and can be used for merging purposes etc.
    # the gtfanno_uid is not used as an index, because the indexing operations used
    # in the following require unique index levels for each annotation, not for each region
    gencode_df, query_df, query_bt = get_input_data(gtf_fp, query_bed, tmpdir_path)
    chromosome_dtype = query_df.Chromosome.dtype

    # %% Produce annotation dataframes for different features, concatenate, sort
    # --------------------------------------------------------------------------

    print("Annotating promoter regions")
    # The annotation dataframes produced in the following steps use the gtfanno_uid as index
    promoter_anno_df = annotate_with_tss(
        gencode_df,
        query_bt,
        distant_cis_regulatory_domain,
        promoter,
        output_cols,
        chromosome_dtype,
        tmpdir_path,
    )

    print("Annotating transcript parts")
    transcript_feature_anno_df = annotate_transcript_parts(
        query_bt, gencode_df, output_cols, chromosome_dtype, tmpdir_path
    )

    print("Merge results")
    # Concatenates two integer indices, leading to misleading duplicates -> reset
    full_annos = pd.concat(
        [promoter_anno_df, transcript_feature_anno_df], axis=0
    ).reset_index(drop=True)
    # Note: not yet sorted! First, we need to make the feature_class categorical

    precedence = pd.Series(
        [
            "Promoter",
            "5'-UTR",
            "3'-UTR",
            "exon",
            "intron",
            # if a feature has no annotation has exon or UTR and it is in a transcript, it is an intron. This if found by annotating query regions against full transcript feature regions and then applying this logic.
            "transcript",
            "DCRD",
            "intergenic",
        ]
    )
    feat_class_cat_dtype = CategoricalDtype(ordered=True, categories=precedence)
    full_annos["feat_class"] = full_annos["feat_class"].astype(feat_class_cat_dtype)

    # This sorting scheme is used for categorizing annotations into primary and secondary
    # As a tie breaker when several annotations for a feature with the same precedence
    # are present, currently the distance is used. Therefore, the distance is included into
    # the sorting scheme. This could be changed into other statistics, and be made available
    # as a parameter
    full_annos_sorted = full_annos.sort_values(
        ["gtfanno_uid", "feat_class", "distance"]
    ).reset_index(drop=True)

    # How many regions can we expect to have multiple annotations?
    # almost all regions have multiple annotations (tested with DMRs against gencode)
    # multi_anno_regions = full_annos_sorted.groupby(['Chromosome', 'Start', 'End']).filter(lambda df: df.shape[0] > 1)

    # %%

    print("Classify annotations")
    full_annos_sorted["feature_rank"] = compute_feature_rank_classification(
        full_annos_sorted
    )
    assert full_annos_sorted['feature_rank'].isin(['primary', 'secondary']).all()
    assert full_annos_sorted['feat_class'].notnull().all()
    # %%

    print("Add intergenic regions")
    # All query regions without annotations are labeled as 'intergenic'
    # Don't add the info columns for regions with annotations, we will fill those
    # with NA via pd.concat
    uids_annotated_regions = full_annos_sorted["gtfanno_uid"].unique()
    intergenic_regions = (
        query_df.loc[~query_df["gtfanno_uid"].isin(uids_annotated_regions), :]
        .copy()
        .assign(
            feat_class=lambda df: pd.Series(
                "intergenic", dtype=feat_class_cat_dtype, index=df.index
            ),
            feature_rank="primary",
        )
    )

    # Merge intergenic regions and regions with annotations
    # the missing annotation columns for the intergenic regions trigger a
    # sort of the columns. Disable with sort=False
    all_regions_annotated = pd.concat(
        [full_annos_sorted, intergenic_regions], sort=False, axis=0
    ).sort_values(["Chromosome", "Start", "End"]).reset_index(drop=True)
    assert all_regions_annotated["feat_class"].dtype == feat_class_cat_dtype
    # Sanity check: when we sort by GRange columns, the gtfanno_uid should be
    # sorted too. All original UIDs should still be present.
    assert (
        all_regions_annotated["gtfanno_uid"].unique() == np.arange(query_df.shape[0])
    ).all()

    print("Save results")
    # TODO: document intron handling, this is also a possible improvement
    primary_annotations = all_regions_annotated.query(
        'feature_rank == "primary"'
    ).copy()
    primary_annotations.loc[
        primary_annotations.feat_class.eq("transcript"), "feat_class"
    ] = "intron"
    assert (
        primary_annotations["gtfanno_uid"].unique() == np.arange(query_df.shape[0])
    ).all()

    # Report some basic stats for further sanity checks
    print("Basic stats for primary annotations")
    primary_anno_freq = (
        primary_annotations["gtfanno_uid"]
        .value_counts()
        .value_counts()
        .to_frame()
        .reset_index()
        .set_axis(["#Primary annotations", "Frequency"], axis=1, inplace=False)
    )
    print(primary_anno_freq)
    feature_freq = (
        primary_annotations["feat_class"]
        .value_counts()
        .to_frame()
        .reset_index()
        .set_axis(["Feature", "Frequency"], axis=1, inplace=False)
    )
    print(feature_freq)

    primary_annotations_bed = trunk_path + "_primary-annotations.bed"
    (
        primary_annotations.rename(columns={"Chromosome": "#Chromosome"}).to_csv(
            primary_annotations_bed, sep="\t", header=True, index=False
        )
    )
    primary_annotations.to_pickle(primary_annotations_bed.replace(".bed", ".p"))

    all_annotations_bed = trunk_path + "_all-annotations.bed"
    (
        all_regions_annotated.rename(columns={"Chromosome": "#Chromosome"}).to_csv(
            all_annotations_bed, sep="\t", header=True, index=False
        )
    )
    all_regions_annotated.to_pickle(all_annotations_bed.replace(".bed", ".p"))


def annotate_transcript_parts(
    query_bt, gencode_df, output_cols, chromosome_dtype, tmpdir_path
) -> pd.DataFrame:
    """

    Args:
        query_bt: BedTool of query dataframe
        gencode_df: gencode_df, as provided by get_input_data
        output_cols: all output columns in final order
        tmpdir_path: path to TemporaryDirectory (results created in this folder
           are not cleaned up by this function directly)
        chromosome_dtype: CategoricalDtype for the chromosome. The categories
            should be strings.

    Returns:
        Dataframe with all query region to transcript feature annotations found

    Notes:
        - this function assumes that the gene annotation sets gene
        boundaries based on transcript boundaries. Therefore, the gene
        feature is not queried (because it is already covered by the
        transcript feature). If gene boundaries may lie outside of the
        transcript boundaries, the gene feature should be added as a
        separate feature to be looked up in this function.
        - the following feature are handled and are therefore expected to be available:
          - transcript, exon, UTR
        - start_codon, stop_codon and other features are currently not
        handled. Currently, transcript parts are only annotated if the
        query region center lies within the transcript feature. If features
        such as start_codons would be added, the code would need to be
        changed to not require the query region center to overlap these special features.
        - A region which overlaps the end of a gene body, but not with its center,
          is currently left unannotated. This is debatable and could be changed.
    """

    # Introns are implicitely annotated, if a region is annotated to a transcript,
    # but not to a exon or UTR
    transcript_features = ("transcript", "exon", "UTR")

    # Annotate UTRs as 5' or 3'
    # --------------------------------------------------------------------------
    transcript_features_df = gencode_df.loc[
        gencode_df["feature"].isin(transcript_features), :
    ].reset_index(drop=True)
    transcript_features_df = expand_gtf_attributes(transcript_features_df)
    is_plus_transcript = transcript_features_df.feat_strand.eq("+")
    is_utr_or_transcript = transcript_features_df.feature.isin(["UTR", "exon"])
    plus_strand_utrs_exon_df = transcript_features_df.loc[
        is_plus_transcript & is_utr_or_transcript, :
    ].sort_values(["transcript_id", "Start", "End"], ascending=True)
    minus_strand_utrs_exons_df = transcript_features_df.loc[
        ~is_plus_transcript & is_utr_or_transcript, :
    ].sort_values(["transcript_id", "Start", "End"], ascending=False)

    print("UTR classification")
    plus_features_with_utr_distinction_ser = _assign_utr_location(
        plus_strand_utrs_exon_df
    )
    minus_features_with_utr_distinction_ser = _assign_utr_location(
        minus_strand_utrs_exons_df
    )
    transcript_features_df.loc[
        plus_features_with_utr_distinction_ser.index, "feature"
    ] = plus_features_with_utr_distinction_ser
    transcript_features_df.loc[
        minus_features_with_utr_distinction_ser.index, "feature"
    ] = minus_features_with_utr_distinction_ser

    # Intersect with query regions
    # --------------------------------------------------------------------------
    transcript_features_fp = tmpdir_path.joinpath("transcript_parts.gtf")
    # bedtools seems to loose the first additional column for a gtf format file?
    # to be safe, add the attribute annotations again at the end
    transcript_features_df.iloc[:, 0:9].to_csv(
        transcript_features_fp, sep="\t", header=False, index=False
    )
    transcript_features_bt = BedTool(str(transcript_features_fp))
    transcript_features_intersect_bt = query_bt.intersect(
        transcript_features_bt, wa=True, wb=True
    )
    col_names = ["Chromosome", "Start", "End", "gtfanno_uid"] + [
        "feat_chrom",
        "source",
        "feat_class",
        "feat_start",
        "feat_end",
        "score",
        "feat_strand",
        "frame",
        "attribute",
    ]
    transcript_feature_annos_df = pd.read_csv(
        transcript_features_intersect_bt.fn,
        sep="\t",
        header=None,
        names=col_names,
        dtype={"Chromosome": chromosome_dtype},
    )

    # Add additional stats
    # --------------------------------------------------------------------------
    feat_size = transcript_feature_annos_df.eval("feat_end - feat_start")
    overlap_start = transcript_feature_annos_df.eval("Start - feat_start").where(
        lambda ser: ser.gt(0), 0
    )
    overlap_end = transcript_feature_annos_df.eval("End - feat_start").where(
        lambda ser: ser.lt(feat_size), feat_size
    )
    overlap_size = overlap_end - overlap_start
    region_size = transcript_feature_annos_df.eval("End - Start")
    # region center
    transcript_feature_annos_df["center"] = transcript_feature_annos_df.eval(
        "Start + (End - Start + 1) / 2 - 1"
    )
    transcript_feature_annos_df["feat_center"] = transcript_feature_annos_df.eval(
        "feat_start + (feat_end - feat_start + 1) / 2 - 1"
    )

    transcript_feature_annos_df.loc[
        lambda df: df["feat_strand"] == "+", "distance"
    ] = transcript_feature_annos_df.eval("feat_center - center")
    transcript_feature_annos_df.loc[
        lambda df: df["feat_strand"] == "-", "distance"
    ] = transcript_feature_annos_df.eval("center - feat_center")

    # TODO: is this correct at the boundaries of the interval?
    transcript_feature_annos_df["has_center"] = np.ceil(
        transcript_feature_annos_df["distance"].abs()
    ).le(np.ceil(feat_size / 2))
    transcript_feature_annos_df["perc_feature"] = overlap_size / feat_size
    transcript_feature_annos_df["perc_region"] = overlap_size / region_size
    transcript_feature_annos_df = expand_gtf_attributes(transcript_feature_annos_df)
    transcript_feature_annos_df = transcript_feature_annos_df[output_cols]

    # Restrict to center overlaps
    # --------------------------------------------------------------------------
    # This part will need to be handled with more precision when features such
    # as start_codons are added, where center-overlap should not be required
    # Also, this disards regions at the end of the gene body which overlap the body
    # but not with their center, which may be suboptimal
    transcript_feature_annos_df = transcript_feature_annos_df.loc[
        transcript_feature_annos_df.has_center, :
    ]
    return transcript_feature_annos_df


def get_input_data(
    gtf_fp, query_bed, tmpdir_path
) -> Tuple[pd.DataFrame, pd.DataFrame, BedTool]:
    """Provide gencode as df, query as df and BedTool

    For the query table, only the Coordinate columns are used. As a fourth
    column, a running uid for each region is added.
    """

    strand_dtype = CategoricalDtype(["+", "-"], ordered=True)

    gencode_df = pd.read_csv(
        gtf_fp,
        sep="\t",
        header=None,
        comment="#",
        names=[
            "feat_chrom",
            "source",
            "feature",
            "Start",
            "End",
            "score",
            "feat_strand",
            "frame",
            "attribute",
        ],
        dtype={
            "feat_chrom": str,
            "Start": "i8",
            "End": "i8",
            "feat_strand": strand_dtype,
        },
    )
    assert gencode_df.eval("Start <= End").all()
    assert (
        gencode_df["feat_strand"].isin(["+", "-"]).all()
    ), "distance computations require that the strand is defined for all features"

    query_df = pd.read_csv(
        query_bed,
        comment="#",
        names=["Chromosome", "Start", "End"],
        usecols=[0, 1, 2],
        header=None,
        sep="\t",
        dtype={"Chromosome": str},
    )
    query_df["gtfanno_uid"] = np.arange(query_df.shape[0])
    query_df["Chromosome"] = pd.Categorical(query_df["Chromosome"], ordered=True)
    assert query_df.eval("Start <= End").all()

    query_bed_with_region_ids = tmpdir_path / "query_bed_with_region_ids.bed"
    query_df.to_csv(query_bed_with_region_ids, header=False, index=False, sep="\t")
    query_bt = BedTool(str(query_bed_with_region_ids))

    return gencode_df, query_df, query_bt


# %%
def annotate_with_tss(
    gencode_df,
    query_bt,
    distant_cis_regulatory_domain,
    promoter,
    output_cols,
    chromosome_dtype,
    tmpdir_path,
) -> pd.DataFrame:
    """Compute TSS annotations

    Find all queries whose center is within the promoter or the DCRD.

    The intervals specified by the TSS distances for the promoter and DCRD
    are closed. If a promoter and DCRD boundary overlap, the region is annotated
    as DCRD.
    """

    # Implementation notes
    # This may be easier to achieve with pybedtools.featurefuncs.TSS or bedtools slop
    # bedtools slop does not provide straightforward reduction of transcripts to TSS, though
    # pybedtools.featurefuncs.TSS has a bug

    # First, we define an area around the TSS where we want to look for overlap
    # with query regions. We look for queries where the center of the query lies within the area
    # covered by the promoter or dcrd
    # So looking for regions which have any overlap with the outer boundaries
    # of this area will cache all candidates!

    # The upper DCRD boundary should be smaller than the upper promoter boundary
    # This is already asserted at the beginning of the anno function, but let's do it again
    # the upper bound is the boundary which is more upstream of the feature
    # both boundaries may be upstream or downstream of the feature
    assert distant_cis_regulatory_domain[0] < promoter[0]
    upper_bound = distant_cis_regulatory_domain[0]
    lower_bound = max(promoter[1], distant_cis_regulatory_domain[1])

    # TSS slop operation to create the TSS area defined above
    transcripts = gencode_df.query('feature == "transcript"').copy()
    transcripts["TSS"] = -1
    transcripts["feat_start"] = -1
    transcripts["feat_end"] = -1
    transcripts["feat_class"] = "TSS_area"
    transcripts = expand_gtf_attributes(transcripts)
    on_plus_strand = transcripts["feat_strand"] == "+"
    transcripts.loc[on_plus_strand, "TSS"] = transcripts.loc[on_plus_strand, "Start"]
    transcripts.loc[~on_plus_strand, "TSS"] = transcripts.loc[~on_plus_strand, "End"]
    transcripts.loc[on_plus_strand, "feat_start"] = (
        transcripts.loc[on_plus_strand, "TSS"] + upper_bound
    )
    transcripts.loc[on_plus_strand, "feat_end"] = (
        transcripts.loc[on_plus_strand, "TSS"] + lower_bound
    )
    # Note that for minus strand features, feat_start and feat_end are still
    # specified as plus strand coordinates
    transcripts.loc[~on_plus_strand, "feat_start"] = (
        transcripts.loc[~on_plus_strand, "TSS"] - lower_bound
    )
    transcripts.loc[~on_plus_strand, "feat_end"] = (
        transcripts.loc[~on_plus_strand, "TSS"] - upper_bound
    )
    transcripts = transcripts.sort_values(
        ["feat_chrom", "feat_start", "feat_end", "TSS"]
    )
    transcripts.loc[transcripts["feat_start"].lt(0), "feat_start"] = 0
    transcripts_cols = [
        "feat_chrom",
        "feat_start",
        "feat_end",
        "TSS",
        "feat_strand",
        "feat_class",
        "gene_name",
        "gene_id",
        "transcript_id",
        "appris_principal_score",
    ]

    # Intersect the TSS areas with the regions_bed.
    # - tss which are further than the downstream and upstream bp away from a region are ignored
    #   and will not be reported
    # - multiple hits per region are possible
    # TODO: check that this works
    all_tss_area_bed = tmpdir_path.joinpath("tss-area-slop.bed")
    transcripts[transcripts_cols].to_csv(
        all_tss_area_bed, sep="\t", header=False, index=False
    )
    all_tss_area_bt = BedTool(str(all_tss_area_bed))
    tss_intersect_bt = query_bt.intersect(all_tss_area_bt, wa=True, wb=True)
    tss_anno = pd.read_csv(
        tss_intersect_bt.fn,
        sep="\t",
        names=["Chromosome", "Start", "End", "gtfanno_uid"] + transcripts_cols,
        dtype={"Chromosome": chromosome_dtype},
    )

    # Add the remaining required output columns and bring in correct order
    tss_anno["perc_feature"] = np.nan
    tss_anno["perc_region"] = np.nan
    # Add distance of query center to TSS
    # has_center is set to False, evaluation whether a region is a promoter
    # or DCRD region will be done in subsequent step
    tss_anno["center"] = tss_anno.eval("Start + (End - Start)/2")
    tss_anno["feat_center"] = np.nan
    tss_anno["has_center"] = False

    tss_anno["distance"] = np.nan
    tss_anno.loc[tss_anno["feat_strand"] == "+", "distance"] = tss_anno.eval(
        "center - TSS"
    )
    tss_anno.loc[tss_anno["feat_strand"] == "-", "distance"] = -tss_anno.eval(
        "center - TSS"
    )
    assert tss_anno["distance"].notnull().all()
    tss_anno = tss_anno[output_cols]

    # classify into proximal and distal cis regulatory regions
    # discard other tss area overlaps
    is_proximal_promoter = tss_anno.eval(f"{promoter[0]} <= distance <= {promoter[1]}")
    is_distant_cis_regulatory_domain = tss_anno.eval(
        f"{distant_cis_regulatory_domain[0]} "
        f"<= distance"
        f" <= {distant_cis_regulatory_domain[1]}"
    )
    tss_anno["feat_class"] = np.nan
    tss_anno.loc[is_proximal_promoter, "feat_class"] = "Promoter"
    tss_anno.loc[is_proximal_promoter, "has_center"] = True
    tss_anno.loc[is_distant_cis_regulatory_domain, "feat_class"] = "DCRD"
    tss_anno.loc[is_distant_cis_regulatory_domain, "has_center"] = True
    # Discard regions without annotation
    # These are regions which overlap with TSS area, but not with their center
    tss_anno = tss_anno.loc[~tss_anno["feat_class"].isna(), :]

    return tss_anno


# %%


def expand_gtf_attributes(df):
    """Extracts metadata from attributes and saves in new dataframe columns
    The following fields are parsed and provided as new columns:
    - gene_id
    - transcript_id
    - gene_name
    - appris_principal_score (float). Features without a principal score are
      assigned a score of 0
    """
    df["gene_id"] = df["attribute"].str.extract('gene_id "(.*?)";')
    df["transcript_id"] = df["attribute"].str.extract('transcript_id "(.*?)";')
    df["gene_name"] = df["attribute"].str.extract('gene_name "(.*?)";')
    df["appris_principal_score"] = (
        df["attribute"].str.extract('tag "appris_principal_(\d)";').astype(float)
    )
    df["appris_principal_score"] = df["appris_principal_score"].fillna(0)
    return df


def _assign_utr_location(utrs_exon_df: pd.DataFrame) -> pd.Series:
    """

    Args:
        utrs_exon_df: Start End feature transcript_id
            sorted on ["transcript_id", "Start", "End"], ascending or descending depending on strand
            df can only contain one strand at a time
            must only contain exons and UTRs in features column

    Returns:
        features series with 5' and 3' UTR distinction

    """

    assert (np.unique(utrs_exon_df["feature"]) == np.array(["UTR", "exon"])).all()

    features_with_utr_5_and_3_prime_distinction = pd.Series(
        "exon", index=utrs_exon_df.index, dtype=object
    )
    # for debugging
    # print_counter = 0
    # from IPython.display import display

    for unused_transcript_id, transcript_group_df in utrs_exon_df.groupby(
        "transcript_id"
    ):
        transcript_interval_is_present_twice = transcript_group_df.duplicated(
            subset=["Start", "End"], keep=False
        )
        non_utr_exon_found = (
            False  # exon is not completely a UTR, it may be partially a UTR
        )
        three_prime_utr_found = False
        # for debugging
        # five_prime_utr_found = True
        for row_idx in range(transcript_group_df.shape[0]):
            curr_feature_name = transcript_group_df.iloc[row_idx].loc["feature"]
            if curr_feature_name == "exon":
                if transcript_interval_is_present_twice.iat[row_idx]:
                    # the entire exon is a UTR
                    # there will be a duplicate entry for the UTR feature, sorted either to occur right before or right after this entry
                    continue
                else:
                    # at least a part of this exon is not a UTR
                    non_utr_exon_found = True
                    # defensive check: after having called a 3' UTR, no more non UTR exons may occur
                    if three_prime_utr_found:
                        raise ValueError("Misclassification of UTR occured")
            elif curr_feature_name == "UTR":
                if non_utr_exon_found:
                    features_with_utr_5_and_3_prime_distinction.at[
                        transcript_group_df.index[row_idx]
                    ] = "3'-UTR"
                    # Used for defensive check above
                    three_prime_utr_found = True
                else:
                    features_with_utr_5_and_3_prime_distinction.at[
                        transcript_group_df.index[row_idx]
                    ] = "5'-UTR"
                    # for printing during debugging
                    # five_prime_utr_found = True
            else:
                raise ValueError(
                    "Unexpected feature encountered during UTR classification"
                )
        # for debugging, print transcripts containing UTRs
        # if (three_prime_utr_found or five_prime_utr_found) and transcript_group_df['feat_strand'].iat[0] == '-':
        #     display(transcript_group_df[["Start", "End", "feature", "feat_strand"]])
        #     display(
        #         features_with_utr_5_and_3_prime_distinction.loc[
        #             transcript_group_df.index
        #         ]
        #     )
        #     print_counter += 1
        #     if print_counter > 50:
        #         raise ValueError()

    return features_with_utr_5_and_3_prime_distinction


def _classify_annos(feat_class_ser, gene_names_ser):
    """Called from groupby-transform: classify annotations into primary and secondary for one region

    Algorithm:
    - the sorting scheme is leveraged to find the highest priority candidate:
      the highest priority feature will be at the beginning of the series,
      and within the feature, the lowest distance will be at the beginning
    """

    # Only one annotation for the region -> must be primary
    if feat_class_ser.shape[0] == 1:
        rank_ser = pd.Series("primary", index=feat_class_ser.index)
        return rank_ser

    rank_ser = pd.Series("secondary", index=feat_class_ser.index)
    top_class = feat_class_ser.iat[0]
    is_top_class = feat_class_ser.eq(top_class)

    # Only one annotation for the top class -> make the top annotation primary
    if is_top_class.sum() == 1:
        rank_ser.loc[is_top_class] = "primary"
        return rank_ser

    # If we get to here, we have more than one top feature hit
    # Do we have top feature hits for several genes, e.g. are we in
    # the promoter of two different genes (perhaps one on the plus
    # and one on the minus strand)?
    gene_anno = gene_names_ser.loc[feat_class_ser.index]
    if gene_anno.nunique() == 1:
        # We have only one gene
        # The data are presorted by the distance, so we can just mark the top hit
        # Note that this may e.g. be promoters from several
        # transcripts for the same gene, and we don't prefer the
        # primary transcript here... This could easily be
        # implemented by sorting by the appris_principal score
        rank_ser.iat[0] = "primary"
        return rank_ser

    # If we get to here, we have multiple genes for the top feature
    # For each gene, take the best top feature hit - currently the
    # one with the smallest distance. So currently no preference for
    # the principal transcript.
    def tag_first_element(ser):
        ser = ser.copy()
        ser.iat[0] = "primary"
        return ser

    rank_ser.loc[is_top_class] = (
        rank_ser.loc[is_top_class].groupby(gene_anno).transform(tag_first_element)
    )

    return rank_ser


def compute_feature_rank_classification(annos_sorted):
    """Compute feature rank classification ('primary' vs 'secondary')

    Args:
        annos_sorted:
            - The dataframe MUST BE SORTED on ["gtfanno_uid", "feat_class", "distance"]
            - annotation dataframe following gtfanno conventions, used columns:
                - gtfanno_uid
                - feat_class
                - distance
                - gene_name
            - feature class must be ordered categorical with most important feature as first element, *following the same precedence as used elsewhere in the annotation workflow*
            - in the main workflow, this is called on the full annotation dataframe with all annotations (except for 'intergenic', which is an implementation detail). Note that this function also works with a subset of the full annotations!

    Returns:
        rank classification series with values 'primary' or 'secondary'
        in the original dataframe order
    """

    assert annos_sorted['feat_class'].dtype.name == 'category'

    # assert that input dataframe is appropriately sorted
    pd.testing.assert_frame_equal(
        annos_sorted,
        annos_sorted.sort_values(["gtfanno_uid", "feat_class", "distance"]),
    )
    # and has a unique index (just a defensive measure to assert that we can check that the index order did not change in the result series)
    assert not annos_sorted.index.has_duplicates

    res = annos_sorted.groupby(["gtfanno_uid"], sort=True)["feat_class"].transform(
        _classify_annos, gene_names_ser=annos_sorted["gene_name"]
    )
    # due to gtfanno_uid sorting in groupby + transform, the element order should remain unchanged compared to the original dataframe

    pd.testing.assert_index_equal(res.index, annos_sorted.index)

    return res
