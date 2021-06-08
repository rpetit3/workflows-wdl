version 1.0

task mercury_json {
    input {
        File metadata
    }

    command <<<
        python3 << CODE
        import json
        with open('~{metadata}', 'rt') as fh:
        out = json.load(fh)
        for k,v in out.items():
            with open(k, 'wt') as k_fh:
                k_fh.write(f'{v}\n')
        CODE
    >>>

    output {
        String sample_name            = read_string("sample_name")
        String platform               = read_string("platform")
        String submission_id          = read_string("submission_id")
        String collection_date        = read_string("collection_date")
        String gisaid_submitter       = read_string("gisaid_submitter")
        String Authors                = read_string("Authors")
        String bioproject_accession   = read_string("bioproject_accession")
        String iso_state              = read_string("iso_state")
        String iso_country            = read_string("VERSION")
        String iso_continent          = read_string("iso_continent")
        String collecting_lab         = read_string("collecting_lab")
        String collecting_lab_address = read_string("collecting_lab_address")
        String submitting_lab         = read_string("submitting_lab")
        String subLab_address         = read_string("subLab_address")
        String gender                 = read_string("gender")
        String patient_age            = read_string("patient_age")
  }

    runtime {
        docker: "python:slim"
        memory: "1 GB"
        cpu: 1
        disks: "local-disk 20 HDD"
        dx_instance_type: "mem1_ssd1_v2_x2"
    }
}