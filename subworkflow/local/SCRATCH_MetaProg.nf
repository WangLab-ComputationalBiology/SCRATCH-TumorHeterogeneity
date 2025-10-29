#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include {
  METAPROG_LEIDEN;
  METAPROG_NMF_PREP;
  METAPROG_NMF;
  METAPROG_POST
} from '../../modules/local/MetaProg-DE/main.nf'   

workflow SCRATCH_MetaProg {

  take:
    ch_input_seurat   // value channel from parent

  main:
    // QMD paths defined in nextflow.config
    def nb_leiden = Channel.fromPath(params.leiden_qmd,    checkIfExists: true)
    def nb_prep   = Channel.fromPath(params.nmfprep_qmd,   checkIfExists: true)
    def nb_nmf    = Channel.fromPath(params.nmf_qmd,       checkIfExists: true)
    def nb_post   = Channel.fromPath(params.postnmf_qmd,   checkIfExists: true)
    
        // 1) Leiden (single)
    def (le_report, le_figs, le_data, le_rds) =
      METAPROG_LEIDEN( ch_input_seurat.combine(nb_leiden) )

    // ---- ADD these 2 lines (explicit barrier) ----
    def LEI_TICK = le_rds.map{ 1 }.take(1)          // turn file into a single tick
    def gated_seurat = ch_input_seurat.combine(LEI_TICK).map { seurat, _ -> seurat }
    // ----------------------------------------------

    // 2) Preprocessing (single) — starts only after Leiden finishes, uses original Seurat
    // def (prep_report, prep_figs, prep_data, prep_preprocessed) =
    //   METAPROG_NMF_PREP( gated_seurat.combine(nb_prep) )
    def (prep_report, prep_figs, prep_data, prep_preprocessed) =
      METAPROG_NMF_PREP( gated_seurat.combine(nb_prep) )

    // ---- Fan-out: one job per *_preprocessed.rds ----
    // Split list-of-paths -> single paths, then attach a sample_id
    // fan-out: build (sample_id, file, notebook) tuples
    def nmf_pairs = prep_preprocessed.flatten().map { p ->
      def sid = p.baseName.replaceFirst(/_preprocessed\.rds$/, '')
      tuple(sid, p)
    }
    def (nmf_report, nmf_figs, nmf_data, nmf_rds) =
      METAPROG_NMF( nmf_pairs.combine(nb_nmf) )


    // 4) PostNMF (single) — waits for all NMF fits
    // METAPROG_POST( nmf_rds.collect().combine(nb_post) )
    // METAPROG_POST( nb_post.combine( nmf_rds.collect() ) )
    def post_in = nmf_rds.collect().map { fits_list ->
      tuple( file(params.postnmf_qmd), fits_list )
    }
    METAPROG_POST( post_in )

}
