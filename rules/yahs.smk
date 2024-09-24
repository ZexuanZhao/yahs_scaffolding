rule yahs:
    conda:
        os.path.join(workflow.basedir,"envs/yahs.yaml")
    input:
        assembly = os.path.join(out_dir, "assembly", "assembly_before_scaffolding.fasta"),
        bam = os.path.join(out_dir,"hic_mapping", "combined.filtered.addRG.purged.sorted.bam")
    output:
        assembly = os.path.join(out_dir, "assembly", "assembly_after_scaffolding_scaffolds_final.fa"),
        bin = os.path.join(out_dir, "assembly", "assembly_after_scaffolding.bin"),
        agp = os.path.join(out_dir, "assembly", "assembly_after_scaffolding_scaffolds_final.agp")
    params:
        prefix = os.path.join(out_dir, "assembly", "assembly_after_scaffolding")
    log:
        os.path.join(out_dir, "log", "yahs.log")
    threads:
        threads
    shell:
        """
        yahs \
            --no-contig-ec \
            --no-mem-check \
            -o {params.prefix} \
            {input.assembly} {input.bam} \
            2>{log}
        """

rule juicer_pre:
    conda:
        os.path.join(workflow.basedir,"envs/yahs.yaml")
    input:
        bin = os.path.join(out_dir, "assembly", "assembly_after_scaffolding.bin"),
        agp = os.path.join(out_dir, "assembly", "assembly_after_scaffolding_scaffolds_final.agp"),
        fai = os.path.join(out_dir, "assembly", "assembly_before_scaffolding.fasta.fai")
    output:
        txt = os.path.join(out_dir, "assembly", "assembly_after_scaffolding_scaffolds_final.txt")
    log:
        os.path.join(out_dir, "log", "create_files_for_manual_curation_1.log")
    params:
        prefix = os.path.join(out_dir, "assembly", "assembly_after_scaffolding_scaffolds_final")
    threads:
        threads
    shell:
        """
        juicer pre -a -o {params.prefix} {input.bin} {input.agp} {input.fai} >{log} 2>&1
        """

rule contact_matrix:
    conda:
        os.path.join(workflow.basedir,"envs/yahs.yaml")
    input:
        txt = os.path.join(out_dir, "assembly", "assembly_after_scaffolding_scaffolds_final.txt"),
        log = os.path.join(out_dir, "log", "create_files_for_manual_curation_1.log")
    output:
        os.path.join(out_dir, "assembly", "assembly_after_scaffolding_scaffolds_final.hic")
    log:
        os.path.join(out_dir, "log", "create_files_for_manual_curation_2.log")
    params:
        temp_file = os.path.join(out_dir, "assembly", "hic.part"),
        mem = config["mem"],
        juicer_tools_jar = config["juicer_tools_jar"]
    shell:
        """
        (java -jar -Xmx{params.mem} {params.juicer_tools_jar} pre {input.txt} {params.temp_file} 1>{log} 2>{log} <(cat {input.log} | grep PRE_C_SIZE | awk '{{print $2" "$3}}')) && (mv {params.temp_file} {output})
        """
