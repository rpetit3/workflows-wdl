version 1.0

# Workflows from Theiagen's public_health_viral_genomics
# Source: https://github.com/theiagen/public_health_viral_genomics
import "https://raw.githubusercontent.com/theiagen/public_health_viral_genomics/v1.4.3/workflows/wf_mercury_batch.wdl" as mercury_batch
import "https://raw.githubusercontent.com/theiagen/public_health_viral_genomics/v1.4.3/workflows/wf_mercury_se_prep.wdl" as mercury_se_prep
import "https://raw.githubusercontent.com/theiagen/public_health_viral_genomics/v1.4.3/workflows/wf_titan_clearlabs.wdl" as clearlabs
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
        call clearlabs.titan_clearlabs as titan { 
            input:
                samplename = samplename,
                clear_lab_fastq = r1
        }
    }

    call mercury_json_prep.mercury_json {
        input:
            metadata = metadata
    }

    if (mercury_json.samplename == samplename) {
        call mercury_se_prep.mercury_se_prep {   
            input: 
                samplename = samplename,
                reads      = r1,

                # From Titan
                sequence               = titan.assembly_fasta,
                seq_platform           = titan.seq_platform,
                assembly_method        = titan.assembly_method,

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
            vadr_num_alerts                 = titan.vadr_num_alerts
    }

    output {
        # Titan outputs
        String? seq_platform               = titan.seq_platform
        File?   dehosted_reads             = titan.dehosted_reads
        Int?    fastqc_raw                 = titan.fastqc_raw
        Int?    fastqc_clean               = titan.fastqc_clean

        String kraken_version              = titan.kraken_version
        Float  kraken_human                = titan.kraken_human
        Float  kraken_sc2                  = titan.kraken_sc2
        String kraken_report               = titan.kraken_report
        Float  kraken_human_dehosted       = titan.kraken_human_dehosted
        Float  kraken_sc2_dehosted         = titan.kraken_sc2_dehosted
        String kraken_report_dehosted      = titan.kraken_report_dehosted

        File   aligned_bam                 = titan.aligned_bam
        File   aligned_bai                 = titan.aligned_bai
        File   variants_from_ref_vcf       = titan.variants_from_ref_vcf
        String artic_version               = titan.artic_version
        File   assembly_fasta              = titan.assembly_fasta
        Int    number_N                    = titan.number_N
        Int    assembly_length_unambiguous = titan.assembly_length_unambiguous
        Int    number_Degenerate           = titan.number_Degenerate
        Int    number_Total                = titan.number_Total
        Float  pool1_percent               = titan.pool1_percent
        Float  pool2_percent               = titan.pool2_percent
        Float  percent_reference_coverage  = titan.percent_reference_coverage
        String assembly_method             = titan.assembly_method

        File   consensus_stats             = titan.consensus_stats
        File   consensus_flagstat          = titan.consensus_flagstat
        Float  meanbaseq_trim              = titan.meanbaseq_trim
        Float  meanmapq_trim               = titan.meanmapq_trim
        Float  assembly_mean_coverage      = titan.assembly_mean_coverage
        String samtools_version            = titan.samtools_version

        String pango_lineage               = titan.pango_lineage
        String pangolin_conflicts          = titan.pangolin_conflicts
        String pangolin_notes              = titan.pangolin_notes
        String pangolin_version            = titan.pangolin_version
        File   pango_lineage_report        = titan.pango_lineage_report
        String pangolin_docker             = titan.pangolin_docker

        File   nextclade_json              = titan.nextclade_json
        File   auspice_json                = titan.auspice_json
        File   nextclade_tsv               = titan.nextclade_tsv
        String nextclade_clade             = titan.nextclade_clade
        String nextclade_aa_subs           = titan.nextclade_aa_subs
        String nextclade_aa_dels           = titan.nextclade_aa_dels
        String nextclade_version           = titan.nextclade_version

        File   vadr_alerts_list            = titan.vadr_alerts_list
        Int    vadr_num_alerts             = titan.vadr_num_alerts
        String vadr_docker                 = titan.vadr_docker

        # Mercury Prep
        File   deID_assembly               = mercury_se_prep.deID_assembly
        File?  genbank_assembly            = mercury_se_prep.genbank_assembly
        File?  genbank_metadata            = mercury_se_prep.genbank_metadata
        File?  gisaid_assembly             = mercury_se_prep.gisaid_assembly
        File?  gisaid_metadata             = mercury_se_prep.gisaid_metadata

        # Mercury Batch
        File?  GenBank_upload_meta         = mercury_batch.GenBank_upload_meta
        File?  GenBank_upload_fasta        = mercury_batch.GenBank_upload_fasta
        File   GenBank_batched_samples     = mercury_batch.GenBank_batched_samples
        File   GenBank_excluded_samples    = mercury_batch.GenBank_excluded_samples
        File?  GISAID_upload_meta          = mercury_batch.GISAID_upload_meta
        File?  GISAID_upload_fasta         = mercury_batch.GISAID_upload_fasta
        File   GISAID_batched_samples      = mercury_batch.GISAID_batched_samples
        File   GISAID_excluded_samples     = mercury_batch.GISAID_excluded_samples
    }
}