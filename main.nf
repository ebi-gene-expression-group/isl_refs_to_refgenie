#!/usr/bin/env nextflow

IRAP_CONFIGS = Channel.fromPath( "${params.irapConfigDir}/*.conf" )

IRAP_CONFIGS
  .filter( ~/.*homo_sapiens.*/ )
  .map{
    r -> tuple(r.baseName.toString().replace('.conf', ''), r)
  }
  .set{
    SPECIES_CONFIG
  }

process derive_assembly {

    input:
        tuple val(species), file(speciesConfig) from SPECIES_CONFIG

    output:
        tuple val(species), stdout, file(speciesConfig) into SPECIES_CONF_WITH_ASSEMBLY

    """
    cat $speciesConfig | grep reference | awk -F'=' '{print \$2}' | awk -F '.' '{print \$2}' | tr -d \'\\n\'
    """
}

process annotate_configline_with_species {

    input:
        tuple val(species), val(assemblyName), file(speciesConfig) from SPECIES_CONF_WITH_ASSEMBLY

    output:
        file("${speciesConfig}.annotated") into ANNOTATED_CONFIGS

    """
    awk -v species=$species -v assembly=$assemblyName '{print species"="assembly"="\$1}' ${speciesConfig} > ${speciesConfig}.annotated
    """
}

ANNOTATED_CONFIGS
    .splitText()
    .map{r -> r.split('\\=') }
    .map{r -> tuple(r[0].toString(), r[1].toString(), r[2].toString(), r[3].toString().trim())}
    .into{
        CONFIGS_FOR_GENOME
        CONFIGS_FOR_CDNA
        CONFIGS_FOR_ANNOTATION
        CONFIGS_FOR_SALMON
    }   

process build_genome {
 
    conda "${baseDir}/envs/refgenie.yml"

    input:
        tuple val(species), val(assembly), val(fileName) from CONFIGS_FOR_GENOME.filter{ it[2] == 'reference'}.map{r -> tuple(r[0], r[1], r[3])}

    output:
        tuple val(species), file(".done") into GENOME_DONE

    """
    filePath=${params.irapDataDir}/reference/$species/$fileName
    build_asset.sh $assembly fasta fasta \$filePath ${params.refgenieDir} 
    """
}

CONFIGS_FOR_ANNOTATION
    .filter{ it[2] == 'gtf_file'}
    .join( GENOME_DONE )
    .map{r -> tuple(r[0], r[1], r[3])}
    .set{
        ANNOTATION_BUILD_INPUTS
    }

process build_annotation {
    
    conda "${baseDir}/envs/refgenie.yml"

    input:
        tuple val(species), val(assembly), val(fileName) from ANNOTATION_BUILD_INPUTS

    output:
        tuple val(species), file(".done") into ANNOTATION_DONE

    """
    ensembl_release=\$(echo $fileName | awk -F'.' '{print \$(NF-2)}' | tr -d \'\\n\') 
    filePath=${params.irapDataDir}/reference/$species/$fileName
    build_asset.sh $assembly ensembl_gtf ensembl_gtf \$filePath ${params.refgenieDir} e\$ensembl_release 
    """
}

process build_cdna {
 
    conda "${baseDir}/envs/refgenie.yml"

    input:
        tuple val(species), val(assembly), val(fileName) from CONFIGS_FOR_CDNA.filter{ it[2] == 'cdna_file'}.map{r -> tuple(r[0], r[1], r[3])}

    output:
        tuple val(species), file(".done") into CDNA_DONE

    """
    ensembl_release=\$(echo $fileName | awk -F'.' '{print \$(NF-2)}' | tr -d \'\\n\') 
    filePath=${params.irapDataDir}/reference/$species/$fileName
    build_asset.sh ${assembly}_cdna fasta fasta \$filePath ${params.refgenieDir} e\$ensembl_release
    """
}

process build_salmon_index {
 
    conda "${baseDir}/envs/refgenie.yml"

    input:
        tuple val(species), val(assembly), val(fileName) from CONFIGS_FOR_SALMON.filter{ it[2] == 'cdna_file'}.map{r -> tuple(r[0], r[1], r[3])}

    output:
        tuple val(species), file(".done") into SALMON_DONE

    """
    ensembl_release=\$(echo $fileName | awk -F'.' '{print \$(NF-2)}' | tr -d \'\\n\') 
    build_asset.sh ${assembly}_cdna salmon_index '' '' ${params.refgenieDir} e\$ensembl_release
    """
}

