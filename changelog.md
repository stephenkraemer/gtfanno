2020-10-17

- add 'intron' and 'intergenic' to the precedence table
- fix: all_regions_annotated['feat_class'] is object dtype (and probably also primary_annos['feat_class']) instead of categorical
  - note that this bug did not affect the correct classification of primary and secondary annotations or the end-result when merging annotations into the one-row-per-feature dataframe
- revert to original precedence order, ie promoter, 5'-UTR, 3'-UTR
- factor out and document compute_feature_rank_classification(annos_sorted)
  - this can also be used with a subset of all annotations, eg to get one annotation
    per region from only transcript part annos (intron/exon/5'-UTR/3'-UTR)

2020-10-16

commit 1a5dcd4

- fix: primary annos df contains secondary annos if one feature is annotated to multiple genes
- fix: distinction between 3' and 5' UTRs is not always correct
- change precedence order (better solution coming up): 5'-UTR before Promoter, to avoid masking 5'-UTRs with promoter annotation
- fix: upstream and downstream interval for Promoter and DCRR was switched on minus strand
- fix: center, distance and has_center calculation for transcript parts is shifted
- fix: bug in distance calculation from TSS
- add some defensive assertions