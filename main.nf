#!/usr/bin/env nextflow

IRAP_CONFIGS = Channel.fromPath( "${params.irapConfigDir}/homo_sapiens.conf" )
SPIKES_GENOME = Channel.fromPath( "${baseDir}/spikes/*/*.fa.gz" ).filter { !it.toString().contains('transcript') }.map{r -> tuple(r.toString().split('/')[-2], r)}
SPIKES_CDNA = Channel.fromPath( "${baseDir}/spikes/*/*.transcripts.fa.gz" ).map{r -> tuple(r.toString().split('/')[-2], r)}
SPIKES_GTF = Channel.fromPath( "${baseDir}/spikes/*/*.gtf.gz" ).filter  { !it.toString().contains('transcript') }.map{r -> tuple(r.toString().split('/')[-2], r)}
SPIKES=SPIKES_GENOME.join(SPIKES_CDNA).join(SPIKES_GTF)

SPIKES_CDNA_GTF = Channel.fromPath( "${baseDir}/spikes/*/*transcripts.gtf.gz" ).map{r -> tuple(r.toString().split('/')[-2], r)}
GENOMES = Channel.fromPath( "${params.islGenomes}" ).splitText().map{r -> r.split()}.filter{ it.size() == 7 }.map{ r -> r*.toString() }
        
ECOLI = Channel.of(['escherichia_coli', "${params.contamination.ecoli.assembly}", file("${params.contamination.ecoli.uri}"), 'default']).first()
FUNGI = Channel.of(['fungi', "${params.contamination.fungi.assembly}", file("${params.contamination.fungi.uri}"), 'default']).first()
VIRUSES = Channel.of(['viruses', "${params.contamination.viruses.assembly}", file("${params.contamination.viruses.uri}"), 'default']).first()

// Get the species from the file- the name is not completely reliable

process find_current_reference_files {

    executor 'local'

    input: 
        file(confFile) from IRAP_CONFIGS

    output:
        tuple file('species.txt'), file('assembly.txt'), file('release.txt'), file('reference.txt'), file('cdna_file.txt'), file('gtf_file.txt'), file('tags.txt') into CURRENT_REF_FILES

    """
    species=\$(grep "species=" $confFile | awk -F'=' '{print \$2}' | tr -d '\\n')
    echo -n "\$species" > species.txt
    fileRoot=${params.irapDataDir}/reference/\$species

    for field in reference cdna_file gtf_file; do
        echo -n \${fileRoot}/\$(grep "\$field=" $confFile | awk -F'=' '{print \$2}' | tr -d '\\n') > \${field}.txt
    done
    detect_current_isl_genome_assembly.sh $confFile ${params.islGenomes} > assembly.txt
    source=\$(grep "^\$(cat species.txt) " "${params.islGenomes}" | awk '{print \$3}') 
    echo -n \$source\$(detect_current_isl_genome_release.sh $confFile) > release.txt

    echo -n 'current' > tags.txt
    """
}

CURRENT_REF_FILES
    .map{ r -> r*.text }
    .filter{ it[0] == 'homo_sapiens' }
    .into{
        CURRENT_REF_FILES_FOR_NEWEST
        CURRENT_REF_FILES_FOR_DOWNSTREAM
    }

process find_newest_reference_files {

    executor 'local'

    input:
        tuple val(species), val(taxId), val(source), val(genomePattern), val(cdnaPattern), val(gtfPattern), val(assembly) from GENOMES.join(CURRENT_REF_FILES_FOR_NEWEST.map{r -> tuple(r[0])})

    output:
        tuple file('species.txt'), file('assembly.txt'), file('release.txt'), file('reference.txt'), file('cdna_file.txt'), file('gtf_file.txt'), file('tags.txt') into NEWEST_REF_FILES
         
    """
    # Just making these files for consistency
    echo -n "$species" > species.txt
    echo -n "$assembly" > assembly.txt

    newestRelease=\$(detect_newest_isl_genome_release.sh $species ${params.irapDataDir} $gtfPattern)
    echo -n "${source}\${newestRelease}" > release.txt

    fileRoot=${params.irapDataDir}/reference/$species
    echo -n \${fileRoot}/\$(basename $genomePattern | sed "s/RELNO/\${newestRelease}/" | sed 's|primary_assembly|toplevel|') > reference.txt
    echo -n \${fileRoot}/\$(basename $gtfPattern | sed "s/RELNO/\${newestRelease}/") > gtf_file.txt

    unversioned_cdna_fasta=\$(basename $cdnaPattern)
    versioned_cdna_fasta=\$(echo -e "\$unversioned_cdna_fasta" | sed "s/.fa.gz/.\${newestRelease}.fa.gz/")
    
    if [ -e "\${fileRoot}/\$versioned_cdna_fasta" ]; then
        echo -n "\${fileRoot}/\$versioned_cdna_fasta" > cdna_file.txt
    else
        echo -n "\${fileRoot}/\$unversioned_cdna_fasta" > cdna_file.txt
    fi

    echo -n 'newest' > tags.txt
    """ 
}

CURRENT_REF_FILES_FOR_DOWNSTREAM
    .concat(NEWEST_REF_FILES.map{ r -> r*.text })
    .map{tuple(it[0], it[1].replace('.', '_'), it[2], file(it[3]), file(it[4]), file(it[5]), it[6])}
    .set{
        REF_FILES
    } 

REF_FILES
    .into{
        REF_FILES_FOR_SPIKES
        REF_FILES_NOT_SPIKES
    }

process add_spikes {

    executor 'local'

    input:
        tuple val(species), val(assembly), val(release), file(referenceFile), file(cdnaFile), file(gtfFile), val(tags), val(spikesName), val(spikesGenome), val(spikesCdna), val(spikesGtf) from REF_FILES_FOR_SPIKES.combine(SPIKES)

    output:
        tuple val(species), val("${assembly}--spikes_${spikesName}"), val(release), file("${species}-${assembly}-${spikesName}.fa.gz"), file("${species}-${assembly}-${release}-${spikesName}.cdna.fa.gz"), file("${species}-${assembly}-${release}-${spikesName}.gtf.gz"), val(tags) into REF_FILES_WITH_SPIKES

    """
    cat $referenceFile $spikesGenome > ${species}-${assembly}-${spikesName}.fa.gz
    cat $cdnaFile $spikesCdna > ${species}-${assembly}-${release}-${spikesName}.cdna.fa.gz
    cat $gtfFile $spikesGtf > ${species}-${assembly}-${release}-${spikesName}.gtf.gz
    """
}

REF_FILES_NOT_SPIKES
    .concat(REF_FILES_WITH_SPIKES)
    .into{
        REF_FILES_FOR_GENOME
        REF_FILES_FOR_POSTGENOME
    }


//// Build the base genome, with and without the genome spikes

// Make a 'genome' containing contamination indices as used by IRAP etc

process make_contamination_fastas {

    input:
        tuple val(ecoliSpecies), val(ecoliAssembly), file("ecoli.fa.gz"), val(tag) from ECOLI
        tuple val(fungiSpecies), val(fungiAssembly), file("fungi.fa.gz"), val(tag) from FUNGI
        tuple val(virusesSpecies), val(virusesAssembly), file("viruses.fa.gz"), val(tag) from VIRUSES

    output:
        tuple val("${ecoliSpecies}-${fungiSpecies}"), val("${ecoliAssembly}-${fungiAssembly}"), file('ecoli_fungi.fa.gz'), val('default') into ECOLI_FUNGI_CONTAMINATION_FASTA
        tuple val("${fungiSpecies}-${virusesSpecies}"), val("${fungiAssembly}-${virusesAssembly}"), file('fungi_viral.fa.gz'), val('default') into FUNGI_VIRUSES_CONTAMINATION_FASTA
        tuple val("${ecoliSpecies}-${virusesSpecies}"), val("${ecoliAssembly}-${virusesAssembly}"), file('ecoli_viral.fa.gz'), val('default') into ECOLI_VIRUSES_CONTAMINATION_FASTA
        tuple val("${ecoliSpecies}-${fungiSpecies}-${virusesSpecies}"), val("${ecoliAssembly}-${fungiAssembly}-${virusesAssembly}"), file('ecoli_fungi_viral.fa.gz'), val('default') into ECOLI_FUNGI_VIRUSES_CONTAMINATION_FASTA

    """
    cat ecoli.fa.gz fungi.fa.gz > ecoli_fungi.fa.gz
    cat fungi.fa.gz viruses.fa.gz > fungi_viral.fa.gz
    cat ecoli.fa.gz viruses.fa.gz > ecoli_viral.fa.gz
    cat ecoli.fa.gz fungi.fa.gz viruses.fa.gz > ecoli_fungi_viral.fa.gz
    """
}

ECOLI_FUNGI_CONTAMINATION_FASTA
    .concat(FUNGI_VIRUSES_CONTAMINATION_FASTA, ECOLI_VIRUSES_CONTAMINATION_FASTA, ECOLI_FUNGI_VIRUSES_CONTAMINATION_FASTA)
    .into{
        CONTAMINATION_GENOMES_FOR_BUILD
        CONTAMINATION_GENOMES_FOR_BOWTIE2
    }

// Group references for the same assembly, e.g. from newest and current
// reference sets so we only build it once. Provided the 'current' tag is
// present this will be set as the default reference for the species.

REF_FILES_FOR_GENOME
    .groupTuple(by: [0,1])
    .map{r -> tuple(r[0], r[1], r[3][0], r[6].unique().join(','))}
    .concat(ECOLI, FUNGI, VIRUSES)
    .concat(CONTAMINATION_GENOMES_FOR_BUILD)
    .map{tuple(it[0], it[1].replace('.', '_'), it[2], it[3])}
    .set{GENOME_BUILD_INPUTS}

process build_genome {
    
    maxForks 1

    memory { 2.GB * task.attempt }

    errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return  task.exitStatus == 130 || task.exitStatus == 137 || task.attempt < 10  ? 'retry': 'ignore' }
    maxRetries 10
 
    conda "${baseDir}/envs/refgenie.yml"

    input:
        tuple val(species), val(assembly), file(filePath), val(additionalTags) from GENOME_BUILD_INPUTS

    output:
        tuple val(species), val(assembly) into GENOME_REFERENCE

    """
    rebuild=''
    if [ "${task.attempt}" -gt 1 ]; then
        rebuild=' -b true'
    fi
    
    # If this is the current, un-spiked geneome for the species, then alias it
    aliases=''
    if [[ $assembly != *"spikes_"* ]] && [[ $additionalTags == *"current"* ]]; then
        aliases=' -l $species'
    fi
    
    build_asset.sh \
        -a ${species}--$assembly \
        -r fasta \
        -f fasta \
        -p $filePath \
        -d ${params.refgenieDir} \
        -t genome\${aliases}\${rebuild} 
    """
}

GENOME_REFERENCE
    .into{
        GENOME_REFERENCE_FOR_POSTGENOME
        GENOME_REFERENCE_FOR_HISAT
        GENOME_REFERENCE_FOR_BOWTIE2
    }

process build_hisat_index {
 
    maxForks 5

    conda "${baseDir}/envs/refgenie.yml"

    memory { 20.GB * task.attempt }

    errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return  task.exitStatus == 130 || task.exitStatus == 137 || task.attempt < 10  ? 'retry': 'ignore' }
    maxRetries 10

    input:
        tuple val(species), val(assembly) from GENOME_REFERENCE_FOR_HISAT

    output:
        tuple val(species), file(".done") into HISAT_DONE

    """
    hisat2_version=\$(cat ${baseDir}/envs/refgenie.yml | grep hisat2 | awk -F'=' '{print \$2}')
    genome_asset="fasta=${species}--${assembly}/fasta:genome"
    
    rebuild=''
    if [ "${task.attempt}" -gt 1 ]; then
        rebuild=' -b true'
    fi

    build_asset.sh \
        -a ${species}--${assembly} \
        -r hisat2_index \
        -d ${params.refgenieDir} \
        -t genome--hisat2_v\${hisat2_version} \
        -s \$genome_asset \${rebuild}
    """
}

// Build bowtie2 indices for contamination genomes specifically

process build_bowtie2_index {
 
    maxForks 5

    conda "${baseDir}/envs/refgenie.yml"

    memory { 20.GB * task.attempt }

    errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return  task.exitStatus == 130 || task.exitStatus == 137 || task.attempt < 10  ? 'retry': 'ignore' }
    maxRetries 10

    input:
        tuple val(species), val(assembly) from GENOME_REFERENCE_FOR_BOWTIE2.join(CONTAMINATION_GENOMES_FOR_BOWTIE2).map{ r -> tuple(r[0], r[1]) }

    output:
        tuple val(species), file(".done") into BOWTIE2_DONE

    """
    bowtie2_version=\$(cat ${baseDir}/envs/refgenie.yml | grep bowtie2 | awk -F'=' '{print \$2}')
    rebuild=''
    if [ "${task.attempt}" -gt 1 ]; then
        rebuild=' -b true'
    fi
    build_asset.sh \
        -a ${species}--${assembly} \
        -r bowtie2_index \
        -d ${params.refgenieDir} \
        -t v\${bowtie2_version}\${rebuild}
    """
}


// The complex cross logic here is just to allow multiple GTFs per assembly
// (which a join I used initially didn't allow).

GENOME_REFERENCE_FOR_POSTGENOME.map{r -> tuple(r[0] + r[1], r).flatten()}
    .cross(REF_FILES_FOR_POSTGENOME.map{r -> tuple(r[0] + r[1], r).flatten()})
    .map{it[1][1..-1]}
    .into{
        GTF_BUILD_INPUTS
        CDNA_BUILD_INPUTS
    }

process build_annotation {
    
    errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return  task.exitStatus == 130 || task.exitStatus == 137 || task.attempt < 10  ? 'retry': 'ignore' }
    maxRetries 10
    
    maxForks 1
    
    conda "${baseDir}/envs/refgenie.yml"

    input:
        tuple val(species), val(assembly), val(version), val(filePath), val(additionalTags) from GTF_BUILD_INPUTS.map{r -> tuple(r[0], r[1], r[2], r[5], r[6])}

    output:
        tuple val(species), file(".done") into ANNOTATION_DONE

    """
    rebuild=''
    if [ "${task.attempt}" -gt 1 ]; then
        rebuild=' -b true'
    fi
    build_asset.sh \
        -a ${species}--${assembly} \
        -r ensembl_gtf \
        -f ensembl_gtf  \
        -p $filePath \
        -d ${params.refgenieDir} \
        -t ${version},${additionalTags}\${rebuild}
    """
}

process build_cdna {
 
    errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return  task.exitStatus == 130 || task.exitStatus == 137 || task.attempt < 10  ? 'retry': 'ignore' }
    maxRetries 10
    
    maxForks 1
    
    conda "${baseDir}/envs/refgenie.yml"

    input:
        tuple val(species), val(assembly), val(version), file(filePath), val(additionalTags) from CDNA_BUILD_INPUTS.map{r -> tuple(r[0], r[1], r[2], r[4], r[6])}

    output:
        tuple val(species), val(assembly), val(version), val(additionalTags) into CDNA_REFERENCE

    """
    rebuild=''
    if [ "${task.attempt}" -gt 1 ]; then
        rebuild=' -b true'
    fi
    build_asset.sh \
        -a ${species}--${assembly} \
        -r fasta_txome \
        -f fasta \
        -p $filePath \
        -d ${params.refgenieDir} \
        -t ${version},${additionalTags} \
        -x cdna_ \${rebuild}
    """
}

CDNA_REFERENCE
    .into{
        CDNA_REFERENCE_FOR_SALMON
        CDNA_REFERENCE_FOR_KALLISTO
    }

process build_salmon_index {
 
    maxForks 5

    conda "${baseDir}/envs/refgenie.yml"

    memory { 20.GB * task.attempt }

    errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return  task.exitStatus == 130 || task.exitStatus == 137 || task.attempt < 10  ? 'retry': 'ignore' }
    maxRetries 10

    input:
        tuple val(species), val(assembly), val(version), val(additionalTags) from CDNA_REFERENCE_FOR_SALMON

    output:
        tuple val(species), file(".done") into SALMON_DONE

    """
    salmon_version=\$(cat ${baseDir}/envs/refgenie.yml | grep salmon | awk -F'=' '{print \$2}')
    cdna_asset="fasta=${species}--${assembly}/fasta_txome:cdna_${version}"
    
    # Append the salmon version to all the input tags
    tags=\$(for at in \$(echo ${version} ${additionalTags} | tr "," "\\n"); do
        echo "\${at}--salmon_v\${salmon_version}"
    done | tr '\\n' ',' | sed 's/,\$//')  
    rebuild=''
    if [ "${task.attempt}" -gt 1 ]; then
        rebuild=' -b true'
    fi
    build_asset.sh \
        -a ${species}--${assembly} \
        -r salmon_index \
        -d ${params.refgenieDir} \
        -t \$tags \
        -x 'cdna_' \
        -s \$cdna_asset \${rebuild}
    """
}

process build_kallisto_index {
 
    maxForks 5

    conda "${baseDir}/envs/refgenie.yml"

    memory { 20.GB * task.attempt }

    errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return  task.exitStatus == 130 || task.exitStatus == 137 || task.attempt < 10  ? 'retry': 'ignore' }
    maxRetries 10

    input:
        tuple val(species), val(assembly), val(version), val(additionalTags) from CDNA_REFERENCE_FOR_KALLISTO

    output:
        tuple val(species), file(".done") into KALLISTO_DONE

    """
    kallisto_version=\$(cat ${baseDir}/envs/refgenie.yml | grep kallisto | awk -F'=' '{print \$2}')
    cdna_asset="fasta=${species}--${assembly}/fasta_txome:cdna_${version}"
    
    # Append the kallisto version to all the input tags
    tags=\$(for at in \$(echo ${additionalTags} | tr "," "\\n"); do
        echo "\${at}--kallisto_v\${kallisto_version}"
    done | tr '\\n' ',' | sed 's/,\$//')  
    
    rebuild=''
    if [ "${task.attempt}" -gt 1 ]; then
        rebuild=' -b true'
    fi
    build_asset.sh \
        -a ${species}--${assembly} \
        -r kallisto_index \
        -d ${params.refgenieDir} \
        -t \$tags \
        -x 'cdna_' \
        -s \$cdna_asset \${rebuild}
    """
}



