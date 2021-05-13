#!/usr/bin/env nextflow

IRAP_CONFIGS = Channel.fromPath( "${params.irapConfigDir}/*.conf" ).map{r -> tuple(r.baseName.toString().replace('.conf', ''), r)}
SPIKES_GENOME = Channel.fromPath( "${baseDir}/spikes/*/*.fa.gz" ).filter { !it.toString().contains('transcript') }.map{r -> tuple(r.toString().split('/')[-2], r)}
SPIKES_CDNA = Channel.fromPath( "${baseDir}/spikes/*/*.transcripts.fa.gz" ).map{r -> tuple(r.toString().split('/')[-2], r)}
SPIKES_GTF = Channel.fromPath( "${baseDir}/spikes/*/*.gtf.gz" ).filter  { !it.toString().contains('transcript') }.map{r -> tuple(r.toString().split('/')[-2], r)}
SPIKES_CDNA_GTF = Channel.fromPath( "${baseDir}/spikes/*/*transcripts.gtf.gz" ).map{r -> tuple(r.toString().split('/')[-2], r)}
GENOMES = Channel.fromPath( "${params.islGenomes}" )

// Extract general info on references from the top-level ISL genome config

GENOMES
  .splitText()
  .filter{ it.contains('homo_sapiens') }
  .map{ r -> r.split() } 
  .map{r -> tuple(r[0].toString(), [ 'reference', 'cdna_file', 'gtf_file' ], r[6].toString(), r[2].toString(),[ r[3].toString(), r[4].toString(), r[5].toString().trim() ])}
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
        tuple val(species), file(speciesConfig) from IRAP_CONFIGS.filter{it[0] == 'homo_sapiens'}

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
    .map{r -> tuple(r[0], r[2], r[4], r[5])}
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
        tuple val(species), val(assembly), val(source), val(pattern) from GENOME_INFO_FOR_CDNA_NEWEST.map{r -> tuple(r[0], r[2], r[3], r[4])}

    output:
        tuple val(species),  val(assembly), file('version.txt'), file('cdna_filename.txt') into CDNA_NEWEST  

    """
    cdna_pattern=\$(basename $pattern | sed 's/RELNO/\\*/' | sed 's/.fa.gz/.\\*.fa.gz/') 
    cdna_fasta=\$(ls ${params.irapDataDir}/reference/$species/\$cdna_pattern | sort -rV | head -n 1)
    
    genome_version=\$(echo -e \$cdna_fasta | awk -F'.' '{print \$(NF-2)}' | tr -d \'\\n\')
    sourcePrefix=\$(echo $source | head -c 1)
    echo -en "\$sourcePrefix\$genome_version" > version.txt     
    
    if [ ! -e \$cdna_fasta ]; then
        echo "\$cdna_fasta file does not exist" 1>&2
        exit 1
    fi    
    echo -n "\$cdna_fasta" > cdna_filename.txt
    """
}

GENOME_INFO_FOR_CDNA_CURRENT
    .join(CURRENT_CONFIGS_FOR_CDNA, by: [0,1])
    .map{r -> tuple(r[0], r[2], r[3], r[4], r[5])}
    .set{
        CURRENT_CDNA_INPUTS
    }

process find_current_cdna {

    executor 'local'

    input:
        tuple val(species), val(assembly), val(source), val(pattern), val(currentFile) from CURRENT_CDNA_INPUTS

    output:
        tuple val(species), val(assembly), file('version.txt'), file('cdna_filename.txt') into CDNA_CURRENT
        
    """
    current_cdna_path=${params.irapDataDir}/reference/${species}/$currentFile
    if [ ! -e "\$current_cdna_path" ]; then
        echo "\$current_cdna_path does not exist" 1>&2
        exit 1
    fi
    echo -n "\$current_cdna_path" > cdna_filename.txt
    
    genome_version=\$(echo -e $currentFile | awk -F'.' '{print \$(NF-2)}' | tr -d \'\\n\')
    sourcePrefix=\$(echo $source | head -c 1)
    echo -en "\$sourcePrefix\$genome_version" > version.txt     
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
        tuple val(species), val(assembly), val(source), val(pattern) from GENOME_INFO_FOR_GTF_NEWEST.map{r -> tuple(r[0], r[2], r[3], r[4])}

    output:
        tuple val(species), val(assembly), file('version.txt'), file('gtf_filename.txt') into GTF_NEWEST  

    """
    gtf_pattern=\$(basename $pattern | sed 's/RELNO/\\*/')
    gtf=\$(ls ${params.irapDataDir}/reference/${species}/\$gtf_pattern | sort -rV | head -n 1)
   
    if [ ! -e \$gtf ]; then
        echo "\$gtf file does not exist" 1>&2
        exit 1
    fi    
    
    echo -n "\$gtf" > gtf_filename.txt
    genome_version=\$(echo -e \$gtf | awk -F'.' '{print \$(NF-2)}' | tr -d \'\\n\')
    sourcePrefix=\$(echo $source | head -c 1)
    echo -en "\$sourcePrefix\$genome_version" > version.txt     
    """
}

GENOME_INFO_FOR_GTF_CURRENT
    .join(CURRENT_CONFIGS_FOR_GTF, by: [0,1])
    .map{r -> tuple(r[0], r[2], r[3], r[4], r[5])}
    .set{
        CURRENT_GTF_INPUTS
    }

process find_current_gtf {

    executor 'local'

    input:
        tuple val(species), val(assembly), val(source), val(pattern), val(currentFile) from CURRENT_GTF_INPUTS

    output:
        tuple val(species), val(assembly), file('version.txt'), file('gtf_filename.txt') into GTF_CURRENT
        
    """
    current_gtf_path=${params.irapDataDir}/reference/${species}/$currentFile
    if [ ! -e "\$current_gtf_path" ]; then
        echo "\$current_gtf_path does not exist" 1>&2
        exit 1
    fi
    echo -n "\$current_gtf_path" > gtf_filename.txt
    
    genome_version=\$(echo -e $currentFile | awk -F'.' '{print \$(NF-2)}' | tr -d \'\\n\')
    sourcePrefix=\$(echo $source | head -c 1)
    echo -en "\$sourcePrefix\$genome_version" > version.txt     
    """
}

// Build the base genome, with and without the genome spikes

REFERENCE_CURRENT
    .map{r -> tuple(r[0], r[1], r[2].text)}
    .into{
        REFERENCE_CURRENT_FOR_BUILD
        REFERENCE_CURRENT_FOR_SPIKES
    }

process add_genome_spikes {
    
    executor 'local'

    input:
        tuple val(species), val(assembly), val(filePath), val(spikesName), file(spikesFile) from REFERENCE_CURRENT_FOR_SPIKES.combine(SPIKES_GENOME)

    output:
        tuple val(species), val("${assembly}-${spikesName}"), file("${assembly}-${spikesName}.fa.gz") into REFERENCE_CURRENT_WITH_SPIKES   
 
    """
    cat $filePath $spikesFile > ${assembly}-${spikesName}.fa.gz
    """
}

REFERENCE_CURRENT_FOR_BUILD
    .concat(REFERENCE_CURRENT_WITH_SPIKES)
    .set{
        GENOME_BUILD_INPUTS
    }

process build_genome {
 
    conda "${baseDir}/envs/refgenie.yml"

    input:
        tuple val(species), val(assembly), val(filePath) from GENOME_BUILD_INPUTS

    output:
        tuple val(species), val(assembly), file(".done") into GENOME_DONE

    """
    build_asset.sh $assembly fasta fasta $filePath ${params.refgenieDir} 
    """
}

// Build GTF for base genome, with and without spikes

GTF_NEWEST
    .concat(GTF_CURRENT)
    .map{r -> tuple(r[0], r[1], r[2].text, r[3].text)}
    .into{
        GTF_FOR_BUILD
        GTF_FOR_SPIKES
    }
    
process add_genome_gtf_spikes {
    
    input:
        tuple val(species), val(assembly), val(version), val(filePath), val(spikesName), file(spikesFile) from GTF_FOR_SPIKES.combine(SPIKES_GTF)

    output:
        tuple val(species), val("${assembly}-${spikesName}"), val(version), file("${assembly}_${version}-${spikesName}.gtf.gz") into GTF_WITH_SPIKES   
 
    """
    cat $filePath $spikesFile > ${assembly}_${version}-${spikesName}.gtf.gz
    """
}

GTF_FOR_BUILD
    .concat(GTF_WITH_SPIKES)
    .join(GENOME_DONE, by: [0,1])
    .map{r -> tuple(r[0], r[1], r[2], r[3])}
    .set{
        GTF_BUILD_INPUTS
    }

process build_annotation {
    
    conda "${baseDir}/envs/refgenie.yml"

    input:
        tuple val(species), val(assembly), val(version), val(filePath) from GTF_BUILD_INPUTS

    output:
        tuple val(species), file(".done") into ANNOTATION_DONE

    """
    build_asset.sh $assembly ensembl_gtf ensembl_gtf $filePath ${params.refgenieDir} $version
    """
}

// Note: for cDNAs we append the version to the genome assembly

CDNA_NEWEST
    .concat(CDNA_CURRENT)
    .map{r -> tuple(r[0], r[1] + '_cdna_' + r[2].text, r[2].text, r[3].text)}
    .into{
        CDNA_FOR_BUILD
        CDNA_FOR_SPIKES
    }

process add_cdna_spikes {
    
    input:
        tuple val(species), val(assembly), val(version), val(filePath), val(spikesName), file(spikesFile) from CDNA_FOR_SPIKES.combine(SPIKES_CDNA)

    output:
        tuple val(species), val("${assembly}-${spikesName}"), val(version), file("${assembly}-${spikesName}.fa.gz") into CDNA_WITH_SPIKES   
 
    """
    echo "assembly is $assembly version is $version" 1>&2
    cat $filePath $spikesFile > ${assembly}-${spikesName}.fa.gz
    """
}

CDNA_FOR_BUILD
    .concat(CDNA_WITH_SPIKES)
    .set{
        CDNA_BUILD_INPUTS
    }

process build_cdna {
 
    conda "${baseDir}/envs/refgenie.yml"

    input:
        tuple val(species), val(assembly), val(version), val(filePath) from CDNA_BUILD_INPUTS

    output:
        tuple val(species), val("${assembly}"), val(version) into CDNA_REFERENCE

    """
    build_asset.sh ${assembly} fasta fasta $filePath ${params.refgenieDir} $version
    """
}

CDNA_REFERENCE
    .into{
        CDNA_REFERENCE_FOR_SALMON
        CDNA_REFERENCE_FOR_KALLISTO
    }

process build_salmon_index {
 
    conda "${baseDir}/envs/refgenie.yml"

    memory { 20.GB * task.attempt }

    errorStrategy { task.exitStatus == 130 || task.exitStatus == 137 || task.attempt < 3  ? 'retry' : 'ignore' }
    maxRetries 10

    input:
        tuple val(species), val(assembly), val(version) from CDNA_REFERENCE_FOR_SALMON

    output:
        tuple val(species), file(".done") into SALMON_DONE

    """
    salmon_version=$(cat ${baseDir}/envs/refgenie.yml | grep salmon | awk -F'=' '{print $2}')
    build_asset.sh ${assembly} salmon_index '' '' ${params.refgenieDir} v\${salmon_version}
    """
}

process build_kallisto_index {
 
    conda "${baseDir}/envs/refgenie.yml"

    memory { 20.GB * task.attempt }

    errorStrategy { task.exitStatus == 130 || task.exitStatus == 137 || task.attempt < 3  ? 'retry' : 'ignore' }
    maxRetries 10

    input:
        tuple val(species), val(assembly), val(version) from CDNA_REFERENCE_FOR_KALLISTO

    output:
        tuple val(species), file(".done") into KALLISTO_DONE

    """
    kallisto_version=$(cat ${baseDir}/envs/refgenie.yml | grep kallisto | awk -F'=' '{print $2}')
    build_asset.sh ${assembly} kallisto_index '' '' ${params.refgenieDir} v\${kallisto_version}
    """
}

