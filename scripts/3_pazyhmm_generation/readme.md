# PAZy HMM Generation
**Author:** Thomas Quigley\
**Date completed:** May 24th, 2026

**Rationale:** The sequence space of the PETadex was generated using PFAM HMMs. The HMMs used were later discovered to not contain all of the catalytic residues required for the core domain of the PETadex enzymes. For analysis of the "completeness" of the catalytic domain of a PETadex sequence, new HMMs were needed.

**Solution:** An HMM was created for each component of the PETadex (42), allowing for the catalytic residues to be more accurately matched. The creation of these PAZy HMMs were analogous to running jackhmmer on the PETadex BlastNR database (sequences containing X's removed), using input MSA's for each component from the PAZy sequences (211 sequences). The difference in this pipeline, is that: for a sequence to be included in the next iteration of the HMM, it must contain the catalytic residues. This ensured that only sequences that are homologous to the catalytic domain were used to refine the HMM.

**Results:** These new PAZy HMMs were successful in returning PETadex sequences that contain amino acids in the catalytic residues, meaning that the sequence is not a fragment of the catalytic domain, but is "complete". The HMMs also return the amino acid identities within each of the catalytic columns.

[Experimental Data Spreadsheet](https://docs.google.com/spreadsheets/d/1xAr6rwhgy4YE14p2JinmmK6SxaCuXtQZr7GovFSJVzM/edit?gid=1033627188#gid=1033627188)

**Scripts used:**
```
├── # pazyhmm_generation ========================================================
├── pazyhmm.build_pazy_hmms-260524.sh        
├── pazyhmm.parse_hmmsearch_output-260524.sh          
├── pazyhmm.transfer_annotations-260524.sh  
```

**Inputs and Outputs:**
```
├── # s3://petadex/logan/pazy_hmms/ ========================================================
├── hmms/        
├── input_msas  
```

## Methods
### 1. Generate PAZy MSA's
**Spreadsheet Sheet IDs:** first_msa_generation, missing_pazy_msa_generation

All PAZy sequences had the predicted signal sequences removed using SignalP5.0. For each component an MSA was generated using muscle 5.3.linux64 align on the PAZy sequences. Sequences that did not contain the catalytic residues were removed from the MSA, and gappy unaligned columns were removed from the N- and C-terminals. HMMs were created from these MSAs (more info below), and any sequence that did not hit to the HMMs was used to create a second MSA per component, annotated as {component}-1.

### 2. Annotating Reference Sequences
**Spreadsheet Sheet IDs:** mechanism_per_gene, reference_structures, reference.csv

The catalytic mechanism and residues were annotated for each "gene" using known literature. These catalytic residues were annotated on reference structures for each component and compilied into a format used with the transfer annotations script.

### 3. HMM Generation
HMMbuild was used to create HMMs from the input PAZy MSAs. HMMalign and the annotation transfer script was used to annotate the match states containing the catalytic residues. The HMMs were searched against the PETadex BlastNR db and the outputs were parsed, resulting in an alignment of PETadex sequences that contain all of the canonical catalytic residues. The resulting alignment included match states from the initial HMM only if that match state had a > 0.5 coverage across the "catalytically complete" sequences in the alignment. The process was repeated a total of three times, with the resulting alignment being used as input into the next cycle.

### 4. Benchmarking
**Spreadsheet Sheet IDs:** c1_hmm_benchmarking, c9_hmm_benchmarking, all_components_benchmarking, blastnr_extracts_no_x_aa

A few different variations of the HMM generation scripts were used:
1. The %coverage of match states to be tranfered to the new HMM.
2. Which residues are deemed as the canonical catalytic residues.
3. Number of HMM generation cycles.

The results showed that there would be minimal change in the sequences returned when the %coverage was changed, although visual inspection of the alignments showd that lower coverages included noise match states, and high coverages resulted in much smaller HMMs. There was not much of a difference when the canonical residue definition was loosened. The number of HMM hits peaked on the 3rd HMM created, while the number of catalytically intact hits peaked at the 4th HMM.

The HMMs were ran on the entire PETadex BlastNR dataset and the number of fragmented vs. complete HMM hits were documented. Any abnormal ratios were investigated.
