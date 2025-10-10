# PETadex index repo

## [PETAdexR](https://github.com/ababaian/petadexR)
`PETadexR` [[link]](https://github.com/ababaian/petadexR) provides `R` based data exploration and analysis tools for plastic‑degrading enzyme sequences. Leveraging petabase‑scale computational biology and high‑throughput functional genomics, PETadexR helps bioinformaticians and biologists identify, visualize, and prioritize novel PETase homologs for experimental validation.

## [PETadex.io](https://github.com/ababaian/petadex.io/) [Under construction]
`PETadex.io` [[link]](https://github.com/ababaian/petadex.io/) is the repository for the `PETadex website.

## `PETadex` Data-lake (s3://petadex/)
`PETadex-logan` sequences are hosted on S3, accessible via `aws cli` or `wget` of URL paths below.
[ Bug Fix: 2025-10-10. Allowed unsigned cli access, added URL access]

```
s3://petadex/logan
├── # PETadex-Logan ========================================================
├── petadex.logan.enz.v1.fa.zst          # PETadex-Logan Enzymes (18.6 GiB)
├── petadex.logan.seq.v1.fa.zst          # Logan Sequence Matches (61.4 GiB)
├── petadex.logan.all_orfs.v1.fa.zst     # Logan search hits (120.3 GiB)
├── # PETadex-nr ===========================================================
├── petadex.nr.enz.fa.zst                # PETadex-nr Enzymes (144 MB)
├── petadex.nr.seq.fa.zst                # NCBI nr hits (207 MB)
└── petadex.nr.seqcluster.tsv.zst        # PETAdex-nr clustering (38 MB)
```

#### PETadex Logan Enzymes (clustered 90%aaid) : `petadex.logan.enz.v1.fa.zst` (18.6 GiB)
Filtered Petadex-Logan sequences (aa) clustered against NCBI `nr` centroids, and then themselves 

- `aws s3 cp s3://petadex/logan/petadex.logan.enz.v1.fa.zst ./` (preferred)
- `wget https://petadex.s3.us-east-1.amazonaws.com/logan/petadex.logan.enz.v1.fa.zst`
- [Download Link](https://petadex.s3.us-east-1.amazonaws.com/logan/petadex.logan.enz.v1.fa.zst)

#### PETadex Logan Sequences: `petadex.logan.seq.v1.fa.zst` (61.4 GiB)
Logan sequences passing HMM-filtering and coverage criteria for inclusion.

- `aws s3 cp s3://petadex/logan/petadex.logan.seq.v1.fa.zst ./` (preferred)
- `wget https://petadex.s3.us-east-1.amazonaws.com/logan/petadex.logan.seq.v1.fa.zst`
- [Download Link](https://petadex.s3.us-east-1.amazonaws.com/logan/petadex.logan.seq.v1.fa.zst)

#### Raw Logan Sequences: `petadex.logan.all_orfs.v1.fa.zst` (120.3 GiB)
Raw `logan` sequences (amino-acid) aligning to `petadex-nr`. For accessing the nucleotide contigs see: [Logan Contigs](https://github.com/IndexThePlanet/Logan/blob/main/Contigs.md)

- `aws s3 cp s3://petadex/logan/petadex.logan.all_orfs.v1.fa.zst ./` (preferred)
- `wget https://petadex.s3.us-east-1.amazonaws.com/logan/petadex.logan.all_orfs.v1.fa.zst`
- [Download Link](https://petadex.s3.us-east-1.amazonaws.com/logan/petadex.logan.all_orfs.v1.fa.zst)

#### PETadex NR Enzymes (clustered 90%aaid) : `petadex.nr.enz.fa.zst` (144 MB)
Filtered and clustered Petadex-nr sequences (aa) from NCBI `nr`. Search query for Logan.

- `aws s3 cp s3://petadex/logan/petadex.nr.enz.fa.zst ./` (preferred)
- `wget https://petadex.s3.us-east-1.amazonaws.com/logan/petadex.nr.enz.fa.zst`
- [Download Link](https://petadex.s3.us-east-1.amazonaws.com/logan/petadex.nr.enz.fa.zst)
