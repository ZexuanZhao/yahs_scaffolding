rule copy_assembly:
    input:
        config["assembly"]
    output:
        os.path.join(out_dir, "assembly", "assembly_before_scaffolding.fasta")
    threads:
        1
    shell:
        """
        cp {input} {output}
        """

rule index_assembly:
    conda:
        os.path.join(workflow.basedir,"envs/arima_hic_mapping.yaml")
    input:
        os.path.join(out_dir,"assembly", "assembly_before_scaffolding.fasta")
    output:
        os.path.join(out_dir,"assembly", "assembly_before_scaffolding.fasta.bwt"),
        os.path.join(out_dir,"assembly", "assembly_before_scaffolding.fasta.fai")
    threads:
        1
    log:
        os.path.join(out_dir, "log", "index.log")
    shell:
        """
        bwa index {input} 2>{log}
        samtools faidx {input} 2>{log}
        """

rule fastp_hic:
    conda:
        os.path.join(workflow.basedir,"envs/arima_hic_mapping.yaml")
    input:
        config["hiC_read1"],
        config["hiC_read2"]
    threads:
        min(threads, 16)
    output:
        expand(os.path.join(out_dir, "trimmed_reads", "{hic_reads_prefix}.clean.fastq.gz"), hic_reads_prefix = hic_reads_prefixs)
    log:
        expand(os.path.join(out_dir, "log", "{hic_reads_prefix}.fastp.log"), hic_reads_prefix = hic_reads_prefixs)
    shell:
        """
        fastp \
            -i {input[0]} -I {input[1]} \
            -o {output[0]} -O {output[1]} \
            --thread {threads} \
            -j /dev/null -h /dev/null \
            2>{log}
        """

rule hic_mapping:
    conda:
        os.path.join(workflow.basedir,"envs/arima_hic_mapping.yaml")
    input:
        bwt = os.path.join(out_dir,"assembly","assembly_before_scaffolding.fasta.bwt"),
        assembly = os.path.join(out_dir,"assembly","assembly_before_scaffolding.fasta"),
        reads = os.path.join(out_dir,"trimmed_reads", "{hic_reads_prefix}.clean.fastq.gz")
    output:
        os.path.join(out_dir,"hic_mapping", "{hic_reads_prefix}.mapped.bam")
    log:
        os.path.join(out_dir, "log", "{hic_reads_prefix}.mapping.log")
    threads:
        threads
    shell:
        """
        bwa mem \
            -t {threads}\
            {input.assembly} \
            {input.reads} 2>{log} | \
            samtools view -@ {threads} -Sb - \
            > {output} 2>{log}
        """

rule filter5end:
    conda:
        os.path.join(workflow.basedir,"envs/arima_hic_mapping.yaml")
    input:
        os.path.join(out_dir,"hic_mapping","{hic_reads_prefix}.mapped.bam")
    output:
        os.path.join(out_dir,"hic_mapping", "{hic_reads_prefix}.mapped.5endFiltered.bam")
    threads:
        5
    log:
        os.path.join(out_dir, "log", "{hic_reads_prefix}.filter5end.log")
    shell:
        """
        samtools view -h -@ {threads} {input} | \
            filter_five_end.pl | \
            samtools view -Sb -@ {threads} - \
            > {output} 2>{log}
        """

rule conbine_and_filter_bams:
    conda:
        os.path.join(workflow.basedir,"envs/arima_hic_mapping.yaml")
    input:
        bams = expand(os.path.join(out_dir,"hic_mapping", "{hic_reads_prefix}.mapped.5endFiltered.bam"),
            hic_reads_prefix=hic_reads_prefixs),
        fai = os.path.join(out_dir,"assembly","assembly_before_scaffolding.fasta.fai")
    output:
         os.path.join(out_dir,"hic_mapping", "combined.filtered.bam")
    threads:
        5
    params:
        mapq_filter=config["mapq_filter"]
    log:
        os.path.join(out_dir, "log", "conbine_and_filter_bams.log")
    shell:
        """
        two_read_bam_combiner.pl {input.bams} samtools {params.mapq_filter} | \
            samtools view -bS -@ {threads} -t {input.fai} - | \
            samtools sort -@ {threads} -o {output} \
            2>{log}
        """

rule add_readGroup:
    conda:
        os.path.join(workflow.basedir, "envs/arima_hic_mapping.yaml")
    input:
        os.path.join(out_dir, "hic_mapping", "combined.filtered.bam")
    output:
        os.path.join(out_dir, "hic_mapping", "combined.filtered.addRG.bam")
    threads:
        1
    log:
        os.path.join(out_dir, "log", "addRG.log")
    shell:
         """
          picard AddOrReplaceReadGroups \
          I={input} \
          O={output} \
          RGID=1 \
          RGLB=lib1 \
          RGPL=ILLUMINA \
          RGPU=unit1 \
          RGSM=HiC \
          2>{log}
         """

rule mark_duplicate:
    conda:
        os.path.join(workflow.basedir,"envs/arima_hic_mapping.yaml")
    input:
        os.path.join(out_dir,"hic_mapping", "combined.filtered.addRG.bam")
    output:
        metric = os.path.join(out_dir,"qc", "markDuplicate.metric.txt"),
        bam = os.path.join(out_dir,"hic_mapping", "combined.filtered.addRG.purged.bam")
    params:
        mem = config["mem"]
    threads:
        2
    log:
        os.path.join(out_dir, "log", "mark_duplicate.log")
    shell:
        """
        picard MarkDuplicates \
            -Xmx{params.mem} -XX:-UseGCOverheadLimit \
            INPUT={input} \
            OUTPUT={output.bam} \
            METRICS_FILE={output.metric} \
            ASSUME_SORTED=TRUE \
            VALIDATION_STRINGENCY=LENIENT\
            REMOVE_DUPLICATES=TRUE \
            2>{log}
        """

rule sort_by_name_bam:
    conda:
        os.path.join(workflow.basedir,"envs/arima_hic_mapping.yaml")
    input:
        os.path.join(out_dir,"hic_mapping", "combined.filtered.addRG.purged.bam")
    output:
         os.path.join(out_dir,"hic_mapping", "combined.filtered.addRG.purged.sorted.bam")
    threads:
        min(threads, 16)
    log:
        os.path.join(out_dir, "log", "sort_by_name_bam.log")
    shell:
        """
        samtools sort -@ {threads} -o {output} -n {input} 2>{log}
        """
