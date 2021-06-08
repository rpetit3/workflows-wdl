version 1.0

# Workflows from Theiagen's public_health_viral_genomics
# Source: https://github.com/theiagen/public_health_viral_genomics
import "https://github.com/theiagen/public_health_viral_genomics/blob/v1.4.3/workflows/wf_mercury_batch.wdl" as mercury_batch
import "https://github.com/theiagen/public_health_viral_genomics/blob/v1.4.3/workflows/wf_mercury_se_prep.wdl" as mercury_se_prep
import "https://github.com/theiagen/public_health_viral_genomics/blob/v1.4.3/workflows/wf_titan_clearlabs.wdl" as titan_clearlabs
import "../tasks/tasks_prep_mercury.wdl" as mercury_json_prep

workflow phvg {
    meta {
        description: "Incorporates the Titan and Mercury workflows into a single run."
        author: "Robert A. Petit III"
        email:  "robert.petit@theiagen.com"
    }

    input {
        String samplename
        String run_id
        String platform
        File   r1
        File?  r2
        File   metadata
    }

    if (platform == "clearlabs") {
        call titan_clearlabs.titan_clearlabs { 
            input:
                samplename = samplename,
                clear_lab_fastq = r1
        }
    }

    call mercury_json_prep.mercury_json {
        input:
            metadata = metadata_json
    }

    if (mercury_json.samplename == samplename) {
        call mercury_se_prep.mercury_se_prep {   
            input: 
                samplename = samplename,
                reads      = r1,

                # From Titan
                sequence               = titan_clearlabs.assembly_fasta,
                seq_platform           = titan_clearlabs.seq_platform,
                assembly_method        = titan_clearlabs.assembly_method,

                # From metadata
                submission_id          = mercury_json.submission_id,
                collection_date        = mercury_json.collection_date,
                gisaid_submitter       = mercury_json.gisaid_submitter,
                iso_state              = mercury_json.iso_state,
                iso_country            = mercury_json.iso_country,
                iso_continent          = mercury_json.iso_continent,
                collecting_lab         = mercury_json.collecting_lab,
                collecting_lab_address = mercury_json.collecting_lab_address,
                bioproject_accession   = mercury_json.bioproject_accession,
                submitting_lab         = mercury_json.submitting_lab,
                subLab_address         = mercury_json.subLab_address,
                Authors                = mercury_json.Authors,
                gender                 = mercury_json.gender,
                patient_age            = mercury_json.patient_age,
        }
    }

    call mercury_batch.mercury_batch { 
        input:
            samplename                      = samplename,
            genbank_single_submission_fasta = mercury_se_prep.genbank_assembly,
            genbank_single_submission_meta  = mercury_se_prep.genbank_metadata,
            gisaid_single_submission_fasta  = mercury_se_prep.gisaid_assembly,
            gisaid_single_submission_meta   = mercury_se_prep.gisaid_metadata,
            vadr_num_alerts                 = titan_clearlabs.vadr.num_alerts
    }

    output {
        File   report                      = align_and_count.report
        File   report_top_hits             = align_and_count.report_top_hits
        String viral_core_version          = align_and_count.viralngs_version

        # Titan_Clearlabs
        String seq_platform                = seq_method
        File   dehosted_reads              = titan_clearlabs.ncbi_scrub_se.read1_dehosted
        Int    fastqc_raw                  = titan_clearlabs.fastqc_se_raw.number_reads
        Int    fastqc_clean                = titan_clearlabs.fastqc_se_clean.number_reads

        String kraken_version              = titan_clearlabs. kraken2_raw.version
        Float  kraken_human                = titan_clearlabs.kraken2_raw.percent_human
        Float  kraken_sc2                  = titan_clearlabs.kraken2_raw.percent_sc2
        String kraken_report               = titan_clearlabs.kraken2_raw.kraken_report
        Float  kraken_human_dehosted       = titan_clearlabs.kraken2_dehosted.percent_human
        Float  kraken_sc2_dehosted         = titan_clearlabs.kraken2_dehosted.percent_sc2
        String kraken_report_dehosted      = titan_clearlabs.kraken2_dehosted.kraken_report

        File   aligned_bam                 = titan_clearlabs.consensus.trim_sorted_bam
        File   aligned_bai                 = titan_clearlabs.consensus.trim_sorted_bai
        File   variants_from_ref_vcf       = titan_clearlabs.consensus.medaka_pass_vcf
        String artic_version               = titan_clearlabs.consensus.artic_pipeline_version
        File   assembly_fasta              = titan_clearlabs.consensus.consensus_seq
        Int    number_N                    = titan_clearlabs.consensus.number_N
        Int    assembly_length_unambiguous = titan_clearlabs.consensus.number_ATCG
        Int    number_Degenerate           = titan_clearlabs.consensus.number_Degenerate
        Int    number_Total                = titan_clearlabs.consensus.number_Total
        Float  pool1_percent               = titan_clearlabs.consensus.pool1_percent
        Float  pool2_percent               = titan_clearlabs.consensus.pool2_percent
        Float  percent_reference_coverage  = titan_clearlabs.consensus.percent_reference_coverage
        String assembly_method             = titan_clearlabs.consensus.artic_pipeline_version

        File   consensus_stats             = titan_clearlabs.stats_n_coverage.stats
        File   consensus_flagstat          = titan_clearlabs.stats_n_coverage.flagstat
        Float  meanbaseq_trim              = titan_clearlabs.stats_n_coverage_primtrim.meanbaseq
        Float  meanmapq_trim               = titan_clearlabs.stats_n_coverage_primtrim.meanmapq
        Float  assembly_mean_coverage      = titan_clearlabs.stats_n_coverage_primtrim.depth
        String samtools_version            = titan_clearlabs.stats_n_coverage.samtools_version

        String pango_lineage               = titan_clearlabs.pangolin2.pangolin_lineage
        String pangolin_conflicts          = titan_clearlabs.pangolin2.pangolin_conflicts
        String pangolin_notes              = titan_clearlabs.pangolin2.pangolin_notes
        String pangolin_version            = titan_clearlabs.pangolin2.version
        File   pango_lineage_report        = titan_clearlabs.pangolin2.pango_lineage_report
        String pangolin_docker             = titan_clearlabs.pangolin2.pangolin_docker

        File   nextclade_json              = titan_clearlabs.nextclade_one_sample.nextclade_json
        File   auspice_json                = titan_clearlabs.nextclade_one_sample.auspice_json
        File   nextclade_tsv               = titan_clearlabs.nextclade_one_sample.nextclade_tsv
        String nextclade_clade             = titan_clearlabs.nextclade_one_sample.nextclade_clade
        String nextclade_aa_subs           = titan_clearlabs.nextclade_one_sample.nextclade_aa_subs
        String nextclade_aa_dels           = titan_clearlabs.nextclade_one_sample.nextclade_aa_dels
        String nextclade_version           = titan_clearlabs.nextclade_one_sample.nextclade_version

        File   vadr_alerts_list            = titan_clearlabs.vadr.alerts_list
        Int    vadr_num_alerts             = titan_clearlabs.vadr.num_alerts
        String vadr_docker                 = titan_clearlabs.vadr.vadr_docker
    }
}