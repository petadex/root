#!/bin/bash
set -e

trap 'echo "Caught signal — killing all child processes..."; kill -- -$$ 2>/dev/null; exit 1' SIGINT SIGTERM

# ───────────────────── #
# BUILD HMMS FROM PAZY  #
# ───────────────────── #

# ─── INPUT PATH CONFIGURATION ─── #
MASTER_PATH="" # fill with master directory
PAZY_MSA_PATH="${MASTER_PATH}/pazy_approach/component_msas/pazy" #input MSA directory
REFERENCE_FASTAA="${MASTER_PATH}/annotated_sequences/individual_fastaas" # directory with representative amino acid sequences
REFERENCE_ANNOTATIONS="${MASTER_PATH}/annotated_sequences/reference.csv" # csv of the annotations to the representative sequence
BLASTNR_PATH="/~petadex_db/final_fastaa/blastnr_components" # directory with the BlastNR fastaas per component
TRANSFER_SCRIPT="${MASTER_PATH}/pazy_approach/scripts/pazyhmm.transfer_annotations-260524.py" # pazyhmm.transfer_annotations-260524.py
PARSE_SCRIPT="${MASTER_PATH}/pazy_approach/scripts/pazyhmm.parse_hmmsearch_output-260524.py" # pazyhmm.parse_hmmsearch_output-260524.py

# ─── WORKING PATH CONFIGURATION ─── #
LOG_PATH="${MASTER_PATH}/logs/pazy_search/components"

# PAZY HMM FILES
PAZY_HMM_PATH="${MASTER_PATH}/component_hmms/pazy"
PAZY_RESULTS_PATH="${MASTER_PATH}/hmm_results/pazy_search"
PAZY_ANALYSIS_PATH="${MASTER_PATH}/catalytic_hits/pazy_search"
PAZY_ANNOTATIONS_PATH="${MASTER_PATH}/annotated_sequences/pazy_hmm_annotations.csv"
PAZY_ALIGNS="${MASTER_PATH}/annotated_sequences/pazy/hmmalign_output"
PAZY_FRAGMENTS="${MASTER_PATH}/pazy_approach/fragments/pazy"

# BLASTNR HMM FILES
NR_MSA_PATH="${MASTER_PATH}/pazy_approach/component_msas/blastnr"
NR_HMM_PATH="${MASTER_PATH}/pazy_approach/component_hmms/blastnr"
NR_RESULTS_PATH="${MASTER_PATH}/pazy_approach/hmm_results/blastnr_search"
NR_ANALYSIS_PATH="${MASTER_PATH}/pazy_approach/catalytic_hits/blastnr_search"
NR_ANNOTATIONS_PATH="${MASTER_PATH}/annotated_sequences/blastnr_hmm_annotations.csv"
NR_ALIGNS="${MASTER_PATH}/annotated_sequences/blastnr/hmmalign_output"
NR_FRAGMENTS="${MASTER_PATH}/pazy_approach/fragments/blastnr"

# FINAL HMM FILES
FINAL_MSA_PATH="${MASTER_PATH}/pazy_approach/component_msas/final"
FINAL_HMM_PATH="${MASTER_PATH}/pazy_approach/component_hmms/final"
FINAL_RESULTS_PATH="${MASTER_PATH}/pazy_approach/hmm_results/final_search"
FINAL_ANALYSIS_PATH="${MASTER_PATH}/pazy_approach/catalytic_hits/final_search"
FINAL_ANNOTATIONS_PATH="${MASTER_PATH}/annotated_sequences/final_hmm_annotations.csv"
FINAL_ALIGNS="${MASTER_PATH}/annotated_sequences/final/hmmalign_output"
FINAL_FRAGMENTS="${MASTER_PATH}/pazy_approach/fragments/final"

# ─── GUARD AGAINST EMPTY VARIABLES ─── #
: \
  ${LOG_PATH:?} \
  ${PAZY_HMM_PATH:?} ${PAZY_RESULTS_PATH:?} ${PAZY_ANALYSIS_PATH:?} \
  ${PAZY_ALIGNS:?} ${PAZY_FRAGMENTS:?} \
  ${NR_MSA_PATH:?} ${NR_HMM_PATH:?} ${NR_RESULTS_PATH:?} ${NR_ANALYSIS_PATH:?} \
  ${NR_ALIGNS:?} ${NR_FRAGMENTS:?} \
  ${FINAL_MSA_PATH:?} ${FINAL_HMM_PATH:?} ${FINAL_RESULTS_PATH:?} ${FINAL_ANALYSIS_PATH:?} \
  ${FINAL_ALIGNS:?} ${FINAL_FRAGMENTS:?}

# ─── REMOVE AND RECREATE WORKING DIRECTORIES ─── #
rm -rf \
    "${LOG_PATH}" \
    "${PAZY_HMM_PATH}" \
    "${PAZY_RESULTS_PATH}" \
    "${PAZY_ANALYSIS_PATH}" \
    "${PAZY_ALIGNS}" \
    "${PAZY_FRAGMENTS}" \
    "${PAZY_ANNOTATIONS_PATH}" \
    "${NR_MSA_PATH}" \
    "${NR_HMM_PATH}" \
    "${NR_RESULTS_PATH}" \
    "${NR_ANALYSIS_PATH}" \
    "${NR_ALIGNS}" \
    "${NR_FRAGMENTS}" \
    "${NR_ANNOTATIONS_PATH}" \
    "${FINAL_MSA_PATH}" \
    "${FINAL_HMM_PATH}" \
    "${FINAL_RESULTS_PATH}" \
    "${FINAL_ANALYSIS_PATH}" \
    "${FINAL_ALIGNS}" \
    "${FINAL_FRAGMENTS}" \
    "${FINAL_ANNOTATIONS_PATH}"

mkdir -p \
    "${LOG_PATH}" \
    "${PAZY_HMM_PATH}" \
    "${PAZY_RESULTS_PATH}" \
    "${PAZY_RESULTS_PATH}/domtbl" \
    "${PAZY_RESULTS_PATH}/aln" \
    "${PAZY_RESULTS_PATH}/aligned_afa" \
    "${PAZY_RESULTS_PATH}/filtered_sto" \
    "${PAZY_RESULTS_PATH}/e_value_csv" \
    "${PAZY_RESULTS_PATH}/fastaa" \
    "${PAZY_ANALYSIS_PATH}" \
    "${PAZY_ALIGNS}" \
    "${PAZY_FRAGMENTS}" \
    "${NR_MSA_PATH}" \
    "${NR_HMM_PATH}" \
    "${NR_RESULTS_PATH}" \
    "${NR_RESULTS_PATH}/domtbl" \
    "${NR_RESULTS_PATH}/aln" \
    "${NR_RESULTS_PATH}/aligned_afa" \
    "${NR_RESULTS_PATH}/filtered_sto" \
    "${NR_RESULTS_PATH}/e_value_csv" \
    "${NR_RESULTS_PATH}/fastaa" \
    "${NR_ANALYSIS_PATH}" \
    "${NR_ALIGNS}" \
    "${NR_FRAGMENTS}" \
    "${FINAL_MSA_PATH}" \
    "${FINAL_HMM_PATH}" \
    "${FINAL_RESULTS_PATH}" \
    "${FINAL_RESULTS_PATH}/domtbl" \
    "${FINAL_RESULTS_PATH}/aln" \
    "${FINAL_RESULTS_PATH}/aligned_afa" \
    "${FINAL_RESULTS_PATH}/filtered_sto" \
    "${FINAL_RESULTS_PATH}/e_value_csv" \
    "${FINAL_RESULTS_PATH}/fastaa" \
    "${FINAL_ANALYSIS_PATH}" \
    "${FINAL_ALIGNS}" \
    "${FINAL_FRAGMENTS}"

echo "# =========================================== #"
echo "$(date) LOOP 1: PAZY build/search/parse + NR HMM build + hhalign"
echo "# =========================================== #"

for alignment in ${PAZY_MSA_PATH}/*.afa; do
    component=$(basename $alignment .afa)
    clean_component=$(echo "$component" | sed 's/-1$//')

    log="${LOG_PATH}/${component}.log"
    (
        echo "# ───────────────────── #"
        echo "$(date) PAZY: Making ${component} HMM:"
        echo "# ───────────────────── #"
        hmmbuild \
            --symfrac 0.9 \
            "${PAZY_HMM_PATH}/${component}.hmm" \
            "${PAZY_MSA_PATH}/${component}.afa"

        echo "# ───────────────────── #"
        echo "$(date) PAZY: Alinging to ${component} reference sequence:"
        echo "# ───────────────────── #"
        hmmalign \
            -o "${PAZY_ALIGNS}/${component}.out" \
            "${PAZY_HMM_PATH}/${component}.hmm" \
            "${REFERENCE_FASTAA}/annotated_${component}.fa"

        echo "# ───────────────────── #"
        echo "$(date) PAZY: Transfering annotations for ${component} reference sequence:"
        echo "# ───────────────────── #"
        python "${TRANSFER_SCRIPT}" \
            "${PAZY_ALIGNS}/${component}.out" \
            "${REFERENCE_ANNOTATIONS}" \
            "${component}" \
            "${PAZY_ANNOTATIONS_PATH}"
            
        echo "# ───────────────────── #"
        echo "$(date) PAZY: Searching ${clean_component} DB:"
        echo "$(date) PAZY: Searching against $(grep -c ">" ${BLASTNR_PATH}/petadex_blastnr_orfs_${clean_component}.fa) sequences"
        echo "# ───────────────────── #"
        hmmsearch \
            --cpu 10 \
            -Z 5569891 \
            -o /dev/null \
            --domtblout ${PAZY_RESULTS_PATH}/domtbl/${component}.domtbl \
            -A ${PAZY_RESULTS_PATH}/aln/${component}.aln \
            --incE 99999 \
            --incdomE 99999 \
            ${PAZY_HMM_PATH}/${component}.hmm \
            ${BLASTNR_PATH}/petadex_blastnr_orfs_${clean_component}.fa

        echo "# ───────────────────── #"
        echo "$(date) PAZY: Parsing ${component} HMMsearch outputs and creating new MSAs:"
        echo "# ───────────────────── #"
        python "${PARSE_SCRIPT}" \
            ${PAZY_RESULTS_PATH} \
            ${component} \
            ${PAZY_ANALYSIS_PATH} \
            ${PAZY_ANNOTATIONS_PATH} \
            ${PAZY_FRAGMENTS}



    ) >> ${log} 2>&1 &

done
wait
echo "# =========================================== #"
echo "$(date) LOOP 1 complete — PAZy HMM used"
echo "# =========================================== #"

echo "# =========================================== #"
echo "$(date) LOOP 2: NR hmmsearch"
echo "# =========================================== #"
for alignment in ${PAZY_MSA_PATH}/*.afa; do
    component=$(basename $alignment .afa)
    log="${LOG_PATH}/${component}.log"
    clean_component=$(echo "$component" | sed 's/-1$//')

    (
        echo "# ───────────────────── #"
        echo "$(date) BlastNR: Making the ${component} BlastNR HMMs:"
        echo "# ───────────────────── #"
        cp "${PAZY_RESULTS_PATH}/filtered_sto/${component}.sto" "${NR_MSA_PATH}/${component}.sto"
        hmmbuild \
            --hand \
            "${NR_HMM_PATH}/${component}.hmm" \
            "${NR_MSA_PATH}/${component}.sto"

        echo "# ───────────────────── #"
        echo "$(date) BlastNR: Alinging to ${component} reference sequence:"
        echo "# ───────────────────── #"
        hmmalign \
            -o "${NR_ALIGNS}/${component}.out" \
            "${NR_HMM_PATH}/${component}.hmm" \
            "${REFERENCE_FASTAA}/annotated_${component}.fa"

        echo "# ───────────────────── #"
        echo "$(date) BlastNR: Transfering annotations for ${component} reference sequence:"
        echo "# ───────────────────── #"
        python "${TRANSFER_SCRIPT}" \
            "${NR_ALIGNS}/${component}.out" \
            "${REFERENCE_ANNOTATIONS}" \
            "${component}" \
            "${NR_ANNOTATIONS_PATH}"
          
        echo "# ───────────────────── #"
        echo "$(date) BlastNR: Searching BlastNR with the ${clean_component} NR HMM:"
        echo "$(date) BlastNR: Searching against $(grep -c ">" ${BLASTNR_PATH}/petadex_blastnr_orfs_${clean_component}.fa) sequences"
        echo "# ───────────────────── #"
        hmmsearch \
            --cpu 10 \
            -Z 5569891 \
            -o /dev/null \
            --domtblout ${NR_RESULTS_PATH}/domtbl/${component}.domtbl \
            -A ${NR_RESULTS_PATH}/aln/${component}.aln \
            --incE 99999 \
            --incdomE 99999 \
            ${NR_HMM_PATH}/${component}.hmm \
            ${BLASTNR_PATH}/petadex_blastnr_orfs_${clean_component}.fa

        echo "# ───────────────────── #"
        echo "$(date) BlastNR: Parsing ${component} HMMsearch with transferred annotations:"
        echo "# ───────────────────── #"
        python "${PARSE_SCRIPT}" \
            ${NR_RESULTS_PATH} \
            ${component} \
            ${NR_ANALYSIS_PATH} \
            ${NR_ANNOTATIONS_PATH} \
            ${NR_FRAGMENTS}

    ) >> ${log} 2>&1 &

done
wait
echo "# =========================================== #"
echo "$(date) LOOP 2 complete — BlastNR HMM used"
echo "# =========================================== #"


echo "# =========================================== #"
echo "$(date) LOOP 3: Final NR hmmsearch"
echo "# =========================================== #"
for alignment in ${PAZY_MSA_PATH}/*.afa; do
    component=$(basename $alignment .afa)
    clean_component=$(echo "$component" | sed 's/-1$//')
    log="${LOG_PATH}/${component}.log"
    (
        echo "# ───────────────────── #"
        echo "$(date) Final: Making the ${component} Final HMMs:"
        echo "# ───────────────────── #"
        cp "${NR_RESULTS_PATH}/filtered_sto/${component}.sto" "${FINAL_MSA_PATH}/${component}.sto"
        hmmbuild \
            --hand \
            "${FINAL_HMM_PATH}/${component}.hmm" \
            "${FINAL_MSA_PATH}/${component}.sto"

        echo "# ───────────────────── #"
        echo "$(date) Final: Alinging to ${component} reference sequence:"
        echo "# ───────────────────── #"
        hmmalign \
            -o "${FINAL_ALIGNS}/${component}.out" \
            "${FINAL_HMM_PATH}/${component}.hmm" \
            "${REFERENCE_FASTAA}/annotated_${component}.fa"

        echo "# ───────────────────── #"
        echo "$(date) Final: Transfering annotations for ${component} reference sequence:"
        echo "# ───────────────────── #"
        python "${TRANSFER_SCRIPT}" \
            "${FINAL_ALIGNS}/${component}.out" \
            "${REFERENCE_ANNOTATIONS}" \
            "${component}" \
            "${FINAL_ANNOTATIONS_PATH}"
          
        echo "# ───────────────────── #"
        echo "$(date) Final: Searching BlastNR with the ${clean_component} NR HMM:"
        echo "$(date) Final: Searching against $(grep -c ">" ${BLASTNR_PATH}/petadex_blastnr_orfs_${clean_component}.fa) sequences"
        echo "# ───────────────────── #"
        hmmsearch \
            --cpu 10 \
            -Z 5569891 \
            -o /dev/null \
            --domtblout ${FINAL_RESULTS_PATH}/domtbl/${component}.domtbl \
            -A ${FINAL_RESULTS_PATH}/aln/${component}.aln \
            --incE 99999 \
            --incdomE 99999 \
            ${FINAL_HMM_PATH}/${component}.hmm \
            ${BLASTNR_PATH}/petadex_blastnr_orfs_${clean_component}.fa

        echo "# ───────────────────── #"
        echo "$(date) Final: Parsing ${component} HMMsearch with transferred annotations:"
        echo "# ───────────────────── #"
        python "${PARSE_SCRIPT}" \
            ${FINAL_RESULTS_PATH} \
            ${component} \
            ${FINAL_ANALYSIS_PATH} \
            ${FINAL_ANNOTATIONS_PATH} \
            ${FINAL_FRAGMENTS}

    ) >> ${log} 2>&1 &

done
wait
echo "# =========================================== #"
echo "$(date) LOOP 3 complete — Final HMM used"
echo "# =========================================== #"
wait
trap - SIGINT SIGTERM  # Clear the trap on normal exit
echo "# =========================================== #"
echo "$(date) All components are finished!"
echo "# =========================================== #"

