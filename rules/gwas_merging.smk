rule merge_association_studies:
    input:
        assoc = expand("/association_studies/{sample}.assoc", sample=sample_tab.sample_name),
    output:
        combined = "results/aggregated_gwas_data.tsv",
        report = "results/gwas_report.html"
    log:
            "logs/merge_assoc.log"
    threads: 10
    conda:  "../wrappers/merge_gwas/env.yaml"
    script: "../wrappers/merge_gwas/script.py"