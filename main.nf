#!/usr/bin/env nextflow

IRAP_CONFIGS = Channel.fromPath( "${params.irapConfigDir}/*.conf" )
SPIKES_GENOME = Channel.fromPath( "${baseDir}/spikes/*/*.fa.gz" ).filter { !it.toString().contains('transcript') }.map{r -> tuple(r.toString().split('/')[-2], r)}
SPIKES_CDNA = Channel.fromPath( "${baseDir}/spikes/*/*.transcripts.fa.gz" ).map{r -> tuple(r.toString().split('/')[-2], r)}
SPIKES_GTF = Channel.fromPath( "${baseDir}/spikes/*/*.gtf.gz" ).filter  { !it.toString().contains('transcript') }.map{r -> tuple(r.toString().split('/')[-2], r)}
SPIKES_CDNA_GTF = Channel.fromPath( "${baseDir}/spikes/*/*transcripts.gtf.gz" ).map{r -> tuple(r.toString().split('/')[-2], r)}
GENOMES = Channel.fromPath( "${params.islGenomes}" ).splitText().map{r -> r.split()}.filter{ it.size() == 7 }.map{ r -> r*.toString() }
        
ECOLI = Channel.of(['escherichia_coli', "${params.contamination.ecoli.assembly}", file("${params.contamination.ecoli.uri}")]).first()
FUNGI = Channel.of(['fungi', "${params.contamination.fungi.assembly}", file("${params.contamination.fungi.uri}")]).first()
VIRUSES = Channel.of(['viruses', "${params.contamination.viruses.assembly}", file("${params.contamination.viruses.uri}")]).first()

// Get the species from the file- the name is not completely reliable

process find_config_species {

    executor 'local'

    input: 
        file(confFile) from IRAP_CONFIGS

    output:
        tuple stdout, file(confFile) into IRAP_CONFIGS_BY_SPECIES

    """
    grep "species=" $confFile | awk -F'=' '{print \$2}' | tr -d '\\n'
    """
}

// Extract general info on references from the top-level ISL genome config and
// combine with the current config

GENOMES
  .join(IRAP_CONFIGS_BY_SPECIES.unique())
  .into{
    MERGED_CONFIG_FOR_RELEASE
    MERGED_CONFIG_FOR_ANNOTATE
    MERGED_CONFIG_FOR_NEWEST_VERSION
    MERGED_CONFIG_FOR_BUILDS
  }

// Find out what version (Ensembl release etc) we have)

process find_curent_release {

    executor 'local'

    input:
        tuple val(species), file(speciesConfig) from MERGED_CONFIG_FOR_RELEASE.map { r -> tuple(r[0], r[7])}

    output:
        tuple val(species), stdout into CURRENT_VERSIONS

    """
    detect_current_isl_genome_release.sh $speciesConfig
    """
}

process find_newest_release {

    executor 'local'

    input:
        tuple val(species), val(gtfPattern) from MERGED_CONFIG_FOR_NEWEST_VERSION.map{ r -> tuple(r[0], r[5]) }

    output:
        tuple val(species), stdout into NEWEST_VERSIONS

    """
    detect_newest_isl_genome_release.sh $species ${params.irapDataDir} $gtfPattern
    """
}

MERGED_CONFIG_FOR_BUILDS
  .map{r -> tuple(r[0].toString(), [ 'reference', 'cdna_file', 'gtf_file' ], r[6].toString().replaceAll('\\.', '_'), r[2].toString(),[ r[3].toString(), r[4].toString(), r[5].toString().trim() ])}
  .join(CURRENT_VERSIONS)
  .join(NEWEST_VERSIONS)
  .transpose()
  .into{
    GENOME_INFO_FOR_REFERENCE
    GENOME_INFO_FOR_CDNA
    GENOME_INFO_FOR_GTF
  }

// Extract current reference info from the config files
 
process annotate_configline_with_species {
    
    executor 'local'

    input:
        tuple val(species), file(speciesConfig) from MERGED_CONFIG_FOR_ANNOTATE.map { r -> tuple(r[0], r[7])}

    output:
        file("${speciesConfig}.annotated") into ANNOTATED_CONFIGS

    """
    awk -v species=$species '{print species"="\$1}' ${speciesConfig} > ${speciesConfig}.annotated
    """
}

ANNOTATED_CONFIGS
    .splitText()
    .map{r -> r.split('\\=') }
    .filter{ it.length == 3 }
    .map{r -> tuple(r[0].toString(), r[1].toString(), r[2].toString().trim())}
    .into{
        CURRENT_CONFIGS_FOR_GENOME
        CURRENT_CONFIGS_FOR_CDNA
        CURRENT_CONFIGS_FOR_GTF
    }

// Find the current reference

GENOME_INFO_FOR_REFERENCE
    .filter{ it[1] == 'reference'}
    .join(CURRENT_CONFIGS_FOR_GENOME, by: [0,1])
    .map{r -> tuple(r[0], r[2], r[4], r[7])}
    .set{
        CURRENT_REFERENCE_INPUTS
    }

process find_current_reference {

    executor 'local'

    input:
        tuple val(species), val(assembly), val(pattern), val(currentFile) from CURRENT_REFERENCE_INPUTS

    output:
        tuple val(species), val(assembly), file('reference_filename.txt') into REFERENCE_CURRENT
        
    """
    current_reference_path=${params.irapDataDir}/reference/${species}/$currentFile
    if [ ! -e "\$current_reference_path" ]; then
        echo "\$current_reference_path does not exist" 1>&2
        exit 1
    fi
    echo -n "\$current_reference_path" > reference_filename.txt
    """
}

// Find the GTF and CDNA for the newest files ISL has

GENOME_INFO_FOR_CDNA
    .filter{ it[1] == 'cdna_file'}
    .into{
        GENOME_INFO_FOR_CDNA_NEWEST
        GENOME_INFO_FOR_CDNA_CURRENT
    }

process find_newest_cdna {

    executor 'local'

    input:
        tuple val(species), val(configType), val(assembly), val(source), val(pattern), val(currentRelease), val(newestRelease) from GENOME_INFO_FOR_CDNA_NEWEST

    output:
        tuple val(species),  val(assembly), file('version.txt'), file('cdna_filename.txt'), val('newest') optional true into CDNA_NEWEST  

    """
    versioned_cdna_fasta=${params.irapDataDir}/reference/$species/\$(basename $pattern | sed 's/.fa.gz/.${newestRelease}.fa.gz/')
    if [ -e "\$versioned_cdna_fasta" ]; then
        echo "cDNA with explicit version is available"
        echo -n "\$versioned_cdna_fasta" > cdna_filename.txt
        echo -en "${source}${newestRelease}" > version.txt
    fi
    """
}

GENOME_INFO_FOR_CDNA_CURRENT
    .join(CURRENT_CONFIGS_FOR_CDNA, by: [0,1])
    .set{
        CURRENT_CDNA_INPUTS
    }

process find_current_cdna {

    executor 'local'

    input:
        tuple val(species), val(configType), val(assembly), val(source), val(pattern), val(currentRelease), val(newestRelease), val(currentFile) from CURRENT_CDNA_INPUTS

    output:
        tuple val(species), val(assembly), file('version.txt'), file('cdna_filename.txt'), val('current') into CDNA_CURRENT
        
    """
    current_cdna_path=${params.irapDataDir}/reference/${species}/$currentFile
    if [ ! -e "\$current_cdna_path" ]; then
        echo "\$current_cdna_path does not exist" 1>&2
        exit 1
    fi
    echo -n "\$current_cdna_path" > cdna_filename.txt
    echo -en "${source}${currentRelease}" > version.txt     
    """
}

GENOME_INFO_FOR_GTF
    .filter{ it[1] == 'gtf_file'}
    .into{
        GENOME_INFO_FOR_GTF_NEWEST
        GENOME_INFO_FOR_GTF_CURRENT
    }

process find_newest_gtf {

    executor 'local'

    input:
        tuple val(species), val(configType), val(assembly), val(source), val(pattern), val(currentRelease), val(newestRelease) from GENOME_INFO_FOR_GTF_NEWEST

    output:
        tuple val(species), val(assembly), file('version.txt'), file('gtf_filename.txt'), val('newest') into GTF_NEWEST  

    """
    newest_gtf=${params.irapDataDir}/reference/$species/\$(basename $pattern | sed 's/RELNO/${newestRelease}/')
   
    if [ ! -e \$newest_gtf ]; then
        echo "\$newest_gtf file does not exist" 1>&2
        exit 1
    fi    
    
    echo -n "\$newest_gtf" > gtf_filename.txt
    echo -en "${source}${newestRelease}" > version.txt     
    """
}

GENOME_INFO_FOR_GTF_CURRENT
    .join(CURRENT_CONFIGS_FOR_GTF, by: [0,1])
    .set{
        CURRENT_GTF_INPUTS
    }

process find_current_gtf {

    executor 'local'

    input:
        tuple val(species), val(configType), val(assembly), val(source), val(pattern), val(currentRelease), val(newestRelease), val(currentFile) from CURRENT_GTF_INPUTS

    output:
        tuple val(species), val(assembly), file('version.txt'), file('gtf_filename.txt'), val('current') into GTF_CURRENT
        
    """
    current_gtf_path=${params.irapDataDir}/reference/${species}/$currentFile
    if [ ! -e "\$current_gtf_path" ]; then
        echo "\$current_gtf_path does not exist" 1>&2
        exit 1
    fi
    echo -n "\$current_gtf_path" > gtf_filename.txt
    
    echo -en "${source}${currentRelease}" > version.txt     
    """
}

//// Build the base genome, with and without the genome spikes

// Make a 'genome' containing contamination indices as used by IRAP etc

process make_contamination_fastas {

    input:
        tuple val(ecoliSpecies), val(ecoliAssembly), file("ecoli.fa.gz") from ECOLI
        tuple val(fungiSpecies), val(fungiAssembly), file("fungi.fa.gz") from FUNGI
        tuple val(virusesSpecies), val(virusesAssembly), file("viruses.fa.gz") from VIRUSES

    output:
        tuple val("${ecoliSpecies}-${fungiSpecies}"), val("${ecoliAssembly}-${fungiAssembly}"), file('ecoli_fungi.fa.gz') into ECOLI_FUNGI_CONTAMINATION_FASTA
        tuple val("${fungiSpecies}-${virusesSpecies}"), val("${fungiAssembly}-${virusesAssembly}"), file('fungi_viral.fa.gz') into FUNGI_VIRUSES_CONTAMINATION_FASTA
        tuple val("${ecoliSpecies}-${virusesSpecies}"), val("${ecoliAssembly}-${virusesAssembly}"), file('ecoli_viral.fa.gz') into ECOLI_VIRUSES_CONTAMINATION_FASTA
        tuple val("${ecoliSpecies}-${fungiSpecies}-${virusesSpecies}"), val("${ecoliAssembly}-${fungiAssembly}-${virusesAssembly}"), file('ecoli_fungi_viral.fa.gz') into ECOLI_FUNGI_VIRUSES_CONTAMINATION_FASTA

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

REFERENCE_CURRENT
    .map{ r -> tuple( r[0], r[1], file(r[2].text) ) }
    .concat(ECOLI, FUNGI, VIRUSES)
    .concat(CONTAMINATION_GENOMES_FOR_BUILD)
    .map{r -> tuple(r[0], r[0] + '--' + r[1].replace('.', '_'), r[2], 'genome')}
    .into{
        REFERENCE_CURRENT_FOR_BUILD
        REFERENCE_CURRENT_FOR_SPIKES
    }

process add_genome_spikes {
    
    executor 'local'
    
    input:
        tuple val(species), val(assembly), file(filePath), val(additionalTags), val(spikesName), file(spikesFile) from REFERENCE_CURRENT_FOR_SPIKES.combine(SPIKES_GENOME)

    output:
        tuple val(species), val("${assembly}"), file("${assembly}--${spikesName}.fa.gz"), val("genome--spikes_${spikesName}")  into REFERENCE_CURRENT_WITH_SPIKES   
 
    """
    cat $filePath $spikesFile > ${assembly}--${spikesName}.fa.gz
    """
}

process build_genome {
    
    maxForks 1

    memory { 2.GB * task.attempt }

    errorStrategy { task.exitStatus == 130 || task.exitStatus == 137 || task.attempt < 3  ? 'retry' : 'terminate' }
    maxRetries 10
 
    conda "${baseDir}/envs/refgenie.yml"

    input:
        tuple val(species), val(assembly), file(filePath), val(additionalTags) from REFERENCE_CURRENT_FOR_BUILD

    output:
        tuple val(species), val(assembly), val(additionalTags), file(".done") into GENOME_REFERENCE

    """
    build_asset.sh $assembly fasta fasta $filePath ${params.refgenieDir} ${additionalTags} 
    """
}

GENOME_REFERENCE
    .map{r -> tuple(r[0], r[1], r[2])}
    .into{
        GENOME_REFERENCE_FOR_SPIKES_GENOME
        GENOME_REFERENCE_FOR_GTF
        GENOME_REFERENCE_FOR_CDNA
        GENOME_REFERENCE_FOR_HISAT
        GENOME_REFERENCE_FOR_BOWTIE2
    }

// We do the spikes geneome in a separate process so that we can hang it off
// the plain genome

REFERENCE_CURRENT_WITH_SPIKES
    .combine(GENOME_REFERENCE_FOR_SPIKES_GENOME, by: [0,1])
    .map{ r -> r[0..3]}
    .set{
        SPIKES_GENOME_INPUT
    }

process build_genome_with_spikes {
    
    maxForks 1

    memory { 2.GB * task.attempt }

    errorStrategy { task.exitStatus == 130 || task.exitStatus == 137 || task.attempt < 3  ? 'retry' : 'terminate' }
    maxRetries 10
 
    conda "${baseDir}/envs/refgenie.yml"

    input:
        tuple val(species), val(assembly), file(filePath), val(additionalTags) from SPIKES_GENOME_INPUT

    output:
        tuple val(species), val(assembly), val(additionalTags), file(".done") into SPIKES_GENOME_REFERENCE

    """
    build_asset.sh $assembly fasta fasta $filePath ${params.refgenieDir} ${additionalTags} 
    """
}

process build_hisat_index {
 
    maxForks 20

    conda "${baseDir}/envs/refgenie.yml"

    memory { 20.GB * task.attempt }

    errorStrategy { task.exitStatus == 130 || task.exitStatus == 137 || task.attempt < 3  ? 'retry' : 'terminate' }
    maxRetries 10

    input:
        tuple val(species), val(assembly), val(additionalTags) from GENOME_REFERENCE_FOR_HISAT.concat(SPIKES_GENOME_REFERENCE.map{r -> tuple(r[0], r[1], r[2])})

    output:
        tuple val(species), file(".done") into HISAT_DONE

    """
    hisat2_version=\$(cat ${baseDir}/envs/refgenie.yml | grep hisat2 | awk -F'=' '{print \$2}')
    tags=${additionalTags}--hisat\${hisat2_version}
    genome_asset="fasta=${assembly}/fasta:${additionalTags}"
    build_asset.sh ${assembly} hisat2_index '' '' ${params.refgenieDir} ${additionalTags}--hisat2_v\${hisat2_version}
    """
}

// Build bowtie2 indices for contamination genomes specifically

process build_bowtie2_index {
 
    maxForks 20

    conda "${baseDir}/envs/refgenie.yml"

    memory { 20.GB * task.attempt }

    errorStrategy { task.exitStatus == 130 || task.exitStatus == 137 || task.attempt < 3  ? 'retry' : 'terminate' }
    maxRetries 10

    input:
        tuple val(species), val(assembly) from GENOME_REFERENCE_FOR_BOWTIE2.join(CONTAMINATION_GENOMES_FOR_BOWTIE2).map{ r -> tuple(r[0], r[1]) }

    output:
        tuple val(species), file(".done") into BOWTIE2_DONE

    """
    bowtie2_version=\$(cat ${baseDir}/envs/refgenie.yml | grep bowtie2 | awk -F'=' '{print \$2}')
    build_asset.sh ${assembly} bowtie2_index '' '' ${params.refgenieDir} v\${bowtie2_version}
    """
}

// Build GTF for base genome, with and without spikes

GTF_NEWEST
    .concat(GTF_CURRENT)
    .map{r -> tuple(r[0], r[0] + '--' + r[1], r[2].text, file(r[3].text), r[4])}
    .map{r -> tuple(r[0], r[1], r[2], r[3], r[2]+','+r[4]) }
    .into{
        GTF_FOR_BUILD
        GTF_FOR_SPIKES
    }
    
process add_genome_gtf_spikes {
    
    input:
        tuple val(species), val(assembly), val(version), file(filePath), val(additionalTags), val(spikesName), file(spikesFile) from GTF_FOR_SPIKES.combine(SPIKES_GTF)

    output:
        tuple val(species), val(assembly), val(version), file("${assembly}_${version}--${spikesName}.gtf.gz"), stdout into GTF_WITH_SPIKES   
 
    """
    cat $filePath $spikesFile > ${assembly}_${version}--${spikesName}.gtf.gz
    
    # Append the spikes name to all the input tags
    for at in \$(echo ${additionalTags} | tr "," "\\n"); do
        echo "\${at}--spikes_${spikesName}"
    done | tr '\\n' ',' | sed 's/,\$//'   
    """
}

// The complex cross logic here is just to allow multiple GTFs per assembly
// (which a join I used initially didn't allow).

GENOME_REFERENCE_FOR_GTF
    .map{ r -> tuple(r[0] + r[1])}
    .cross( GTF_FOR_BUILD.concat(GTF_WITH_SPIKES).map{ r -> tuple(r[0] + r[1], r[0], r[1], r[2], r[3], r[4]) } )
    .map{ r -> r[1] }
    .map{ r -> tuple( r[1], r[2], r[3], r[4], r[5]) }
    .set{
        GTF_BUILD_INPUTS
    }

process build_annotation {
    
    errorStrategy { task.exitStatus == 130 || task.exitStatus == 137 || task.attempt < 3  ? 'retry' : 'terminate' }
    maxRetries 10
    
    maxForks 1
    
    conda "${baseDir}/envs/refgenie.yml"

    input:
        tuple val(species), val(assembly), val(version), val(filePath), val(additionalTags) from GTF_BUILD_INPUTS

    output:
        tuple val(species), file(".done") into ANNOTATION_DONE

    """
    build_asset.sh $assembly ensembl_gtf ensembl_gtf $filePath ${params.refgenieDir} ${additionalTags}
    """
}

// Note: for cDNAs we append the version to the genome assembly

CDNA_NEWEST
    .concat(CDNA_CURRENT)
    .map{r -> tuple(r[0], r[0] + '--' + r[1], r[2].text, file(r[3].text), r[4])}
    .map{r -> tuple(r[0], r[1], r[2], r[3], r[2]+','+r[4]) }
    .into{
        CDNA_FOR_BUILD
        CDNA_FOR_SPIKES
    }

process add_cdna_spikes {
    
    input:
        tuple val(species), val(assembly), val(version), file(filePath), val(additionalTags), val(spikesName), file(spikesFile) from CDNA_FOR_SPIKES.combine(SPIKES_CDNA)

    output:
        tuple val(species), val(assembly), val(version), file("${assembly}--${spikesName}.fa.gz"), stdout into CDNA_WITH_SPIKES   
 
    """
    cat $filePath $spikesFile > ${assembly}--${spikesName}.fa.gz

    # Append the spikes name to all the input tags
    for at in \$(echo ${additionalTags} | tr "," "\\n"); do
        echo "\${at}--spikes_${spikesName}"
    done | tr '\\n' ',' | sed 's/,\$//'   
    """
}

GENOME_REFERENCE_FOR_CDNA
    .map{ r -> tuple(r[0].toString() + r[1])}
    .cross( CDNA_FOR_BUILD.concat(CDNA_WITH_SPIKES).map{ r -> tuple(r[0] + r[1], r[0], r[1], r[2], r[3], r[4]) } )  
    .map{ r -> r[1] }
    .map{ r -> r[1..5] }
    .set{
        CDNA_BUILD_INPUTS
    }

process build_cdna {
 
    errorStrategy { task.exitStatus == 130 || task.exitStatus == 137 || task.attempt < 3  ? 'retry' : 'terminate' }
    maxRetries 1
    
    maxForks 1
    
    conda "${baseDir}/envs/refgenie.yml"

    input:
        tuple val(species), val(assembly), val(version), file(filePath), val(additionalTags) from CDNA_BUILD_INPUTS

    output:
        tuple val(species), val("${assembly}"), val(version), val(additionalTags) into CDNA_REFERENCE

    """
    build_asset.sh ${assembly} fasta fasta $filePath ${params.refgenieDir} ${additionalTags} cdna_
    """
}

CDNA_REFERENCE
    .into{
        CDNA_REFERENCE_FOR_SALMON
        CDNA_REFERENCE_FOR_KALLISTO
    }

process build_salmon_index {
 
    maxForks 20

    conda "${baseDir}/envs/refgenie.yml"

    memory { 20.GB * task.attempt }

    errorStrategy { task.exitStatus == 130 || task.exitStatus == 137 || task.attempt < 3  ? 'retry' : 'terminate' }
    maxRetries 10

    input:
        tuple val(species), val(assembly), val(version), val(additionalTags) from CDNA_REFERENCE_FOR_SALMON

    output:
        tuple val(species), file(".done") into SALMON_DONE

    """
    salmon_version=\$(cat ${baseDir}/envs/refgenie.yml | grep salmon | awk -F'=' '{print \$2}')
    cdna_asset="fasta=${assembly}/fasta:cdna_\$(echo "${additionalTags}" | awk -F',' '{print \$1}')"
    
    # Append the salmon version to all the input tags
    tags=\$(for at in \$(echo ${additionalTags} | tr "," "\\n"); do
        echo "\${at}--salmon_v\${salmon_version}"
    done | tr '\\n' ',' | sed 's/,\$//')  

    build_asset.sh ${assembly} salmon_index '' '' ${params.refgenieDir} \$tags 'cdna_' \$cdna_asset
    """
}

process build_kallisto_index {
 
    maxForks 20

    conda "${baseDir}/envs/refgenie.yml"

    memory { 20.GB * task.attempt }

    errorStrategy { task.exitStatus == 130 || task.exitStatus == 137 || task.attempt < 3  ? 'retry' : 'terminate' }
    maxRetries 10

    input:
        tuple val(species), val(assembly), val(version), val(additionalTags) from CDNA_REFERENCE_FOR_KALLISTO

    output:
        tuple val(species), file(".done") into KALLISTO_DONE

    """
    kallisto_version=\$(cat ${baseDir}/envs/refgenie.yml | grep kallisto | awk -F'=' '{print \$2}')
    cdna_asset="fasta=${assembly}/fasta:cdna_\$(echo "${additionalTags}" | awk -F',' '{print \$1}')"
    
    # Append the kallisto version to all the input tags
    tags=\$(for at in \$(echo ${additionalTags} | tr "," "\\n"); do
        echo "\${at}--kallisto_v\${kallisto_version}"
    done | tr '\\n' ',' | sed 's/,\$//')  
    
    build_asset.sh ${assembly} kallisto_index '' '' ${params.refgenieDir} \$tags 'cdna_' \$cdna_asset
    """
}

