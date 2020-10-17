2020-10-14

- fix: primary annos df contains secondary annos if one feature is annotated to multiple genes
- fix: distinction between 3' and 5' UTRs is not always correct
- change precedence order (better solution coming up): 5'-UTR before Promoter, to avoid masking 5'-UTRs with promoter annotation
- fix: upstream and downstream interval for Promoter and DCRR was switched on minus strand
- fix: center, distance and has_center calculation for transcript parts is shifted
- fix: bug in distance calculation from TSS
- add some defensive assertions