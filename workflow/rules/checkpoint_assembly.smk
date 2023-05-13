checkpoint check_file_size:
    input:
        os.path.join(config["workdir"], "{PRB}_{EHI}_assembly/", "{EHA}_contigs.fasta")
    output:
        os.path.join(config["workdir"], "{PRB}_{EHI}_assembly/", "assembly_checkpoint")
    shell: 
        """
        test `stat --printf='%s' {input}` -lt 50000000 || touch {output}
        """