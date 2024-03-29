#!/usr/bin/env nextflow

// mode parameter defines if normal and/or splici salmon index is build
mode = params.mode

// IRAP_CONFIGS = Channel.fromPath( "${params.irapConfigDir}/*.conf").filter{it.baseName == 'drosophila_melanogaster'}
IRAP_CONFIGS = Channel.fromPath( "${params.irapConfigDir}/*.conf")
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

    errorStrategy 'finish'
    conda "${baseDir}/envs/refgenie.yml"
    
    input: 
        file(confFile) from IRAP_CONFIGS

    output:
        tuple file('species.txt'), file('assembly.txt'), file('release.txt'), file('reference.txt'), file('cdna_file.txt'), file('gtf_file.txt'), file('tags.txt'), file('build.txt') into CURRENT_REF_FILES

    """
    set +e
    species=\$(grep "species=" $confFile | awk -F'=' '{print \$2}' | tr -d '\\n')
    assembly=\$(detect_current_isl_genome_assembly.sh $confFile ${params.islGenomes})
    release=\$(detect_current_isl_genome_release.sh $confFile "${params.islGenomes}")
    fileRoot=${params.irapDataDir}/reference/\$species
    reference=\$(echo -n \${fileRoot}/\$(grep "^reference=" $confFile | awk -F'=' '{print \$2}' | tr -d '\\n'))    
    cdna_file=\$(echo -n \${fileRoot}/\$(grep "^cdna_file=" $confFile | awk -F'=' '{print \$2}' | tr -d '\\n'))    
    gtf_file=\$(echo -n \${fileRoot}/\$(grep "^gtf_file=" $confFile | awk -F'=' '{print \$2}' | tr -d '\\n'))    
    tag='current'

    env COLUMNS=500 check_refgenie_status.sh "\$species" "\$assembly" "\$release" "\$reference" "\$cdna_file" "\$gtf_file" "\$tag"

    if [ \$? -eq 1 ]; then
        build='true'
    else
        build='false'
    fi
    
    echo -n "\$species" > species.txt
    echo -n "\$assembly" > assembly.txt
    echo -n "\$release" > release.txt
    echo -n "\$reference" > reference.txt
    echo -n "\$cdna_file" > cdna_file.txt
    echo -n "\$gtf_file" > gtf_file.txt
    echo -n "\$tag" > tags.txt
    echo -n "\$build" > build.txt
    """
}

CURRENT_REF_FILES
    .map{ r -> r*.text }
    .into{
        CURRENT_REF_FILES_FOR_NEWEST
        CURRENT_REF_FILES_FOR_DOWNSTREAM
    }

process find_newest_reference_files {

    errorStrategy 'finish'
    conda "${baseDir}/envs/refgenie.yml"

    input:
        tuple val(species), val(taxId), val(source), val(genomePattern), val(cdnaPattern), val(gtfPattern), val(assembly) from GENOMES.join(CURRENT_REF_FILES_FOR_NEWEST.map{r -> tuple(r[0])})

    output:
        tuple file('species.txt'), file('assembly.txt'), file('release.txt'), file('reference.txt'), file('cdna_file.txt'), file('gtf_file.txt'), file('tags.txt'), file('build.txt') into NEWEST_REF_FILES
         
    """
    set +e
    releaseNo=\$(detect_newest_isl_genome_release.sh $species ${params.irapDataDir} $gtfPattern)
    release=${source}\${releaseNo}
    fileRoot=${params.irapDataDir}/reference/$species
    reference=\$(echo -n \${fileRoot}/\$(basename $genomePattern | sed "s/RELNO/\${releaseNo}/" | sed 's|primary_assembly|toplevel|'))
    gtf_file=\$(echo -n \${fileRoot}/\$(basename $gtfPattern | sed "s/RELNO/\${releaseNo}/"))

    # ISL was switched at some point to append the E! release to cDNA files.
    # Ensembl does't do that, but the files do differ between releases. We need
    # to account for versioned and unversioned possibilities.

    unversioned_cdna_fasta=\$(basename $cdnaPattern |  sed "s/RELNO/\${releaseNo}/")
    versioned_cdna_fasta=\$(echo -e "\$unversioned_cdna_fasta" | sed "s/.fa.gz/.\${releaseNo}.fa.gz/")
    
    if [ -e "\${fileRoot}/\$versioned_cdna_fasta" ]; then
        cdna_file=\$(echo -n "\${fileRoot}/\$versioned_cdna_fasta")
    else
        cdna_file=\$(echo -n "\${fileRoot}/\$unversioned_cdna_fasta")
    fi
    tag='newest'

    env COLUMNS=500 check_refgenie_status.sh "$species" "$assembly" "\$release" "\$reference" "\$cdna_file" "\$gtf_file" "\$tag"

    if [ \$? -eq 1 ]; then
        build='true'
    else
        build='false'
    fi

    echo -n "$species" > species.txt
    echo -n "$assembly" > assembly.txt
    echo -n "\$release" > release.txt
    echo -n "\$reference" > reference.txt
    echo -n "\$cdna_file" > cdna_file.txt
    echo -n "\$gtf_file" > gtf_file.txt
    echo -n "\$tag" > tags.txt
    echo -n "\$build" > build.txt
    """ 
}

CURRENT_REF_FILES_FOR_DOWNSTREAM
    .concat(NEWEST_REF_FILES.map{ r -> r*.text })
    .filter{ it[7] == 'true' } 
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
    errorStrategy 'finish'

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

    errorStrategy 'finish'
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
    .into{
        GENOME_ALIAS_INPUTS
        GENOME_BUILD_INPUTS
        SPLICI_BUILD_INPUTS
    }

process build_genome {
    
    memory { 2.GB * task.attempt }

    errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return  task.exitStatus == 130 || task.exitStatus == 137 || task.attempt < 3  ? 'retry': 'ignore' }
    maxRetries 3
 
    conda "${baseDir}/envs/refgenie.yml"

    input:
        tuple val(species), val(assembly), file(filePath), val(additionalTag) from GENOME_BUILD_INPUTS

    output:
        tuple val(species), val(assembly) into GENOME_REFERENCE

    """
    build_asset.sh \
        -a ${species}--$assembly \
        -r fasta \
        -f fasta \
        -p $filePath \
        -d ${params.refgenieDir} \
        -m yes \
        -b true \
        -t genome
    """
}

GENOME_REFERENCE
    .into{
        GENOME_REFERENCE_FOR_COLLECTION
        GENOME_REFERENCE_FOR_POSTGENOME
        GENOME_REFERENCE_FOR_HISAT
        GENOME_REFERENCE_FOR_BOWTIE2
    }

// Build all references in parallel followed by a reduction step before dependent
// processes (e.g. HISAT and Bowtie indexing). This allows us to do lots of
// things in parallel while still having asset dependencies in place at the
// right time

GENOME_REFERENCE_FOR_COLLECTION
    .collect()
    .map { r -> 'collected' }
    .set{
        COLLECTED_REFERENCES
    }

process reduce_genomes {

    conda "${baseDir}/envs/refgenie.yml"
    errorStrategy 'finish'
    
    memory { 20.GB * task.attempt }
    
    maxForks 1
    
    input:
        val(collected) from COLLECTED_REFERENCES

    output:
        val('reduced') into REDUCED_REFERENCES

    """
    refgenie build --reduce -c ${params.refgenieDir}/genome_config.yaml
    """ 
}

// The complex cross logic here is just to allow multiple GTFs per assembly
// (which a join I used initially didn't allow).

GENOME_REFERENCE_FOR_POSTGENOME.map{r -> tuple(r[0] + r[1], r).flatten()}
    .cross(REF_FILES_FOR_POSTGENOME.map{r -> tuple(r[0] + r[1], r).flatten()})
    .map{it[1][1..-1]}
    .into{
        GTF_BUILD_INPUTS
        GTF_SPLICI_BUILD_INPUTS
        CDNA_BUILD_INPUTS
    }

process build_cdna {
 
    errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return  task.exitStatus == 130 || task.exitStatus == 137 || task.attempt < 3  ? 'retry': 'ignore' }
    maxRetries 3
    
    conda "${baseDir}/envs/refgenie.yml"

    input:
        val(reduced) from REDUCED_REFERENCES
        tuple val(species), val(assembly), val(version), file(filePath), val(additionalTags) from CDNA_BUILD_INPUTS.map{r -> tuple(r[0], r[1], r[2], r[4], r[6])}

    output:
        tuple val(species), val(assembly), val(version), val(additionalTags) into CDNA_REFERENCE

    """
    build_asset.sh \
        -a ${species}--${assembly} \
        -r fasta_txome \
        -f fasta \
        -p $filePath \
        -d ${params.refgenieDir} \
        -t ${version},${additionalTags} \
        -m yes \
        -b true \
        -x cdna_ 
    """
}

CDNA_REFERENCE
    .into{
        CDNA_REFERENCE_FOR_COLLECTION
        CDNA_REFERENCE_FOR_SALMON
        CDNA_REFERENCE_FOR_KALLISTO
    }

// Build all cDNAs in parallel followed by a reduction step before dependent
// processes (Salmon and Kallisto indexing) and genome-centric downstream
// processes that would be disrupted by a reduce (hisat and anno). This allows
// us to do lots of things in parallel while still having asset dependencies in
// place at the right time

CDNA_REFERENCE_FOR_COLLECTION
    .collect()
    .map { r -> 'collected' }
    .set{
        COLLECTED_CDNAS
    }

process reduce_cdnas {

    conda "${baseDir}/envs/refgenie.yml"
    errorStrategy 'finish'
    
    memory { 20.GB * task.attempt }
    
    maxForks 1
    
    input:
        val(collected) from COLLECTED_CDNAS

    output:
        val('reduced') into REDUCED_CDNAS

    """
    refgenie build --reduce -c ${params.refgenieDir}/genome_config.yaml
    """ 
}

// build a assets for splici trancritome

process build_splici_txome {
 
    errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return  task.exitStatus == 130 || task.exitStatus == 137 || task.attempt < 3  ? 'retry': 'ignore' }
    maxRetries 3
    
    conda "${baseDir}/envs/refgenie.yml"

    input:
        val(reduced) from REDUCED_REFERENCES
        tuple val(species), val(assembly), file(filePath), val(additionalTag_genome), val(version), file(gtf), val(additionalTags) from SPLICI_BUILD_INPUTS.combine( GTF_SPLICI_BUILD_INPUTS.map{r -> tuple(r[0], r[1], r[2], r[5], r[6])}, by: [0, 1])
       
    output:
        tuple val(species), val(assembly), val(version), val(additionalTags) into SPLICI_REFERENCE
      
    """
    gunzip -f $filePath 
    gunzip -f $gtf
    pyroe make-splici *.fa *.gtf 90 splici --filename-prefix ${species}.${assembly}.
    gzip splici/*.fa
    
    build_asset.sh \
        -a ${species}--${assembly} \
        -r fasta_txome \
        -f fasta \
        -p splici/*.fa.gz \
        -d ${params.refgenieDir} \
        -t ${version},${additionalTags} \
        -m yes \
        -b true \
        -x splici_ 
    """
}

SPLICI_REFERENCE
    .into{
        SPLICI_REFERENCE_FOR_COLLECTION
        SPLICI_REFERENCE_FOR_SALMON
    }
    
    
CDNA_REFERENCE_FOR_SALMON
   .join( SPLICI_REFERENCE_FOR_SALMON)
   .set {
      REFERENCE_FOR_SALMON
      }

SPLICI_REFERENCE_FOR_COLLECTION
    .collect()
    .map { r -> 'collected' }
    .set{
        COLLECTED_SPLICI
    }

process reduce_splici {

    conda "${baseDir}/envs/refgenie.yml"
    errorStrategy 'finish'
    
    memory { 20.GB * task.attempt }
    
    maxForks 1
    
    input:
        val(collected) from COLLECTED_SPLICI

    output:
        val('reduced') into REDUCED_SPLICI

    """
    refgenie build --reduce -c ${params.refgenieDir}/genome_config.yaml
    """ 
}

process build_annotation {
    
    errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return  task.exitStatus == 130 || task.exitStatus == 137 || task.attempt < 3  ? 'retry': 'ignore' }
    maxRetries 3
    
    conda "${baseDir}/envs/refgenie.yml"

    input:
        val(reduced) from REDUCED_REFERENCES
        val(reduced) from REDUCED_CDNAS
        tuple val(species), val(assembly), val(version), val(filePath), val(additionalTags) from GTF_BUILD_INPUTS.map{r -> tuple(r[0], r[1], r[2], r[5], r[6])}

    output:
        tuple val(species), val(assembly), val('none') into ANNOTATION_DONE

    """
    build_asset.sh \
        -a ${species}--${assembly} \
        -r ensembl_gtf \
        -f ensembl_gtf  \
        -p $filePath \
        -d ${params.refgenieDir} \
        -m yes \
        -b true \
        -t ${version},${additionalTags}
    """
}

// In the below, we make everything dependent on the completion of the
// annotation and cDNA builds. This is not because of a functional dependency,
// rather because we can't run 'map' operations alongside 'reduce' operations.
// So we do all the (fairly quick) fasta, gtf and alias building first so we
// can release them for use ASAP, then run the other slow indexing in parallel
// without worrying.

process build_hisat_index {
 
    conda "${baseDir}/envs/refgenie.yml"

    memory { 20.GB * task.attempt }

    errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return  task.exitStatus == 130 || task.exitStatus == 137 || task.attempt < 3  ? 'retry': 'ignore' }
    maxRetries 3

    input:
        val(reduced) from REDUCED_REFERENCES
        val(reduced) from REDUCED_CDNAS
        tuple val(species), val(assembly) from GENOME_REFERENCE_FOR_HISAT

    output:
        tuple val(species), val(assembly), val('none') into HISAT2_DONE

    """
    hisat2_version=\$(cat ${baseDir}/envs/refgenie.yml | grep hisat2 | awk -F'=' '{print \$2}')
    genome_asset="fasta=${species}--${assembly}/fasta:genome"
    
    build_asset.sh \
        -a ${species}--${assembly} \
        -r hisat2_index \
        -d ${params.refgenieDir} \
        -t genome--hisat2_v\${hisat2_version} \
        -m yes \
        -b true \
        -s \$genome_asset
    """
}

// Build bowtie2 indices for contamination genomes specifically

process build_bowtie2_index {
 
    conda "${baseDir}/envs/refgenie.yml"

    memory { 20.GB * task.attempt }

    errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return  task.exitStatus == 130 || task.exitStatus == 137 || task.attempt < 3  ? 'retry': 'ignore' }
    maxRetries 3

    input:
        val(reduced) from REDUCED_REFERENCES
        val(reduced) from REDUCED_CDNAS
        tuple val(species), val(assembly) from GENOME_REFERENCE_FOR_BOWTIE2.join(CONTAMINATION_GENOMES_FOR_BOWTIE2).map{ r -> tuple(r[0], r[1]) }

    output:
        tuple val(species), val(assembly), val('none') into BOWTIE2_DONE

    """
    bowtie2_version=\$(cat ${baseDir}/envs/refgenie.yml | grep bowtie2 | awk -F'=' '{print \$2}')
    build_asset.sh \
        -a ${species}--${assembly} \
        -r bowtie2_index \
        -d ${params.refgenieDir} \
        -m yes \
        -b true \
        -t v\${bowtie2_version}
    """
}

process build_salmon_index {
 
    conda "${baseDir}/envs/refgenie.yml"

    memory { 20.GB * task.attempt }

    errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return  task.exitStatus == 130 || task.exitStatus == 137 || task.attempt < 3  ? 'retry': 'ignore' }
    maxRetries 3

    input:
        val(reduced) from REDUCED_SPLICI
        val(reduced) from REDUCED_CDNAS
        tuple val(species), val(assembly), val(version), val(additionalTags) from REFERENCE_FOR_SALMON

    output:
        tuple val(species), val(assembly), val(version) into SALMON_DONE
        
    script:
    if( mode == 'normal' )
        """
        salmon_version=\$(cat ${baseDir}/envs/refgenie.yml | grep salmon | awk -F'=' '{print \$2}')
        cdna_asset="fasta=${species}--${assembly}/fasta_txome:cdna_${version}"
     
         # Append the salmon version to all the input tags
        tags=\$(for at in \$(echo ${version} ${additionalTags} | tr "," "\\n"); do
                echo "\${at}--salmon_v\${salmon_version}"
        done | tr '\\n' ',' | sed 's/,\$//')  
        build_asset.sh \
         -a ${species}--${assembly} \
         -r salmon_index \
         -d ${params.refgenieDir} \
         -t \$tags \
         -x 'cdna_' \
         -m yes \
         -b true \
         -s \$cdna_asset
        """
    else if( mode == 'splici' )
        """
        salmon_version=\$(cat ${baseDir}/envs/refgenie.yml | grep salmon | awk -F'=' '{print \$2}')
        splici_asset="fasta=${species}--${assembly}/fasta_txome:splici_${version}"
        
        # Append the salmon version to all the input tags
        tags=\$(for at in \$(echo ${version} ${additionalTags} | tr "," "\\n"); do
                 echo "\${at}--salmon_v\${salmon_version}"
        done | tr '\\n' ',' | sed 's/,\$//')  
        build_asset.sh \
        -a ${species}--${assembly} \
        -r salmon_index \
        -d ${params.refgenieDir} \
        -t \$tags \
        -x 'splici_' \
        -m yes \
        -b true \
        -s \$splici_asset
        """
     else if( mode == 'both' )
        """
        #build both the normal and the splici index
        salmon_version=\$(cat ${baseDir}/envs/refgenie.yml | grep salmon | awk -F'=' '{print \$2}')
        cdna_asset="fasta=${species}--${assembly}/fasta_txome:cdna_${version}"
   
        tags=\$(for at in \$(echo ${version} ${additionalTags} | tr "," "\\n"); do
                echo "\${at}--salmon_v\${salmon_version}"
        done | tr '\\n' ',' | sed 's/,\$//')  
        build_asset.sh \
         -a ${species}--${assembly} \
         -r salmon_index \
         -d ${params.refgenieDir} \
         -t \$tags \
         -x 'cdna_' \
         -m yes \
         -b true \
         -s \$cdna_asset


        splici_asset="fasta=${species}--${assembly}/fasta_txome:splici_${version}"
        
        
        tags=\$(for at in \$(echo ${version} ${additionalTags} | tr "," "\\n"); do
                 echo "\${at}--salmon_v\${salmon_version}"
        done | tr '\\n' ',' | sed 's/,\$//')  
        build_asset.sh \
        -a ${species}--${assembly} \
        -r salmon_index \
        -d ${params.refgenieDir} \
        -t \$tags \
        -x 'splici_' \
        -m yes \
        -b true \
        -s \$splici_asset
        """
    else 
        """
        error "Invalid salmon mode: ${mode}"
        """
}

process build_kallisto_index {
 
    conda "${baseDir}/envs/refgenie.yml"

    memory { 20.GB * task.attempt }

    errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return  task.exitStatus == 130 || task.exitStatus == 137 || task.attempt < 3  ? 'retry': 'ignore' }
    maxRetries 3

    input:
        val(reduced) from REDUCED_CDNAS
        tuple val(species), val(assembly), val(version), val(additionalTags) from CDNA_REFERENCE_FOR_KALLISTO

    output:
        tuple val(species), val(assembly), val(version) into KALLISTO_DONE

    """
    kallisto_version=\$(cat ${baseDir}/envs/refgenie.yml | grep kallisto | awk -F'=' '{print \$2}')
    cdna_asset="fasta=${species}--${assembly}/fasta_txome:cdna_${version}"
    
    # Append the kallisto version to all the input tags
    tags=\$(for at in \$(echo ${additionalTags} | tr "," "\\n"); do
        echo "\${at}--kallisto_v\${kallisto_version}"
    done | tr '\\n' ',' | sed 's/,\$//')  
    
    build_asset.sh \
        -a ${species}--${assembly} \
        -r kallisto_index \
        -d ${params.refgenieDir} \
        -t \$tags \
        -x 'cdna_' \
        -m yes \
        -b true \
        -s \$cdna_asset
    """
}

// Collect species-wise outputs and run the refgenie reduce 
// Note: we'd like to be able to run a reduce once we have all the assets for
// each species. But since we can't run reduce operations while maps are
// ongoing, we need to wait until everything is done (hence the collect(). If
// that restriction is removed by the refgeneie bods, we can fix this.

HISAT2_DONE.groupTuple()
    .join(SALMON_DONE.groupTuple())
    .join(KALLISTO_DONE.groupTuple())
    .collect()
    .map { r -> 'collected' }
    .set{
        SPECIES_REDUCTIONS
    }

process reduce {

    conda "${baseDir}/envs/refgenie.yml"
    
    memory { 20.GB * task.attempt }
    
    maxForks 1
    
    input:
        val(collected) from SPECIES_REDUCTIONS

    output:
        val('reduced') into REDUCED_SPECIES

    """
    refgenie build --reduce -c ${params.refgenieDir}/genome_config.yaml
    """ 
}

// Get a list of current aliases, which will be used to determine necessary updates

process get_alias_table {

    conda "${baseDir}/envs/refgenie.yml"
    
    input:
        val('reduced') from REDUCED_SPECIES
   
    output:
        file('alias_table.txt') into ALIAS_TABLE 

    """
    env COLUMNS=500 refgenie alias get > alias_table.txt.tmp
    mv alias_table.txt.tmp alias_table.txt 
    """
}


// By making aliasing dependent on the output of the reduction of indices we
// can use aliasing as a marker of completion for a species

process alias_genomes {

    conda "${baseDir}/envs/refgenie.yml"
    errorStrategy 'ignore'
    
    maxForks 1
    
    input:
        file(aliasTable) from ALIAS_TABLE
        tuple val(species), val(assembly), file(filePath), val(additionalTag) from GENOME_ALIAS_INPUTS

    output:
         tuple val(species), val(assembly), val('none') into ALIAS_DONE

    """
    set +e
    export REFGENIE=${params.refgenieDir}/genome_config.yaml

    # The 'current' assembly for a species will be the defult, so set a
    # species-wide alias for that

    aliases=''
    
    spikePart=''
    echo "$assembly" | grep "spikes" > /dev/null
    if [ \$? -eq -0 ]; then
        assembly=$assembly
        parts=(\${assembly//--/ })
        spikePart=--\${parts[1]} 
    fi
    
    echo -e "$additionalTag" | grep "current" > /dev/null 2>&1
    if [ \$? -eq 0 ]; then
        aliases="$species\${spikePart}"
    fi

    # Non-spiked genomes will be aliased to e.g. homo_sapiens--current. Spiked
    # genomes will be aliased to e.g. homo_sapiens--current--spikes_ercc

    for tag in \$(echo -e "$additionalTag" | sed 's/,/ /g'); do 
        alias="${species}--\$tag"
        if [ -z "\$aliases" ]; then
            aliases=\$alias\${spikePart}
        else
            aliases="\${aliases} \${alias}\${spikePart}"
        fi
    done

    # Now alias this assembly    
    digest=\$(digest_from_alias.sh ${species}--${assembly} $aliasTable)

    for alias in \$aliases; do
        existing_alias_digest=\$(digest_from_alias.sh \$alias ${aliasTable})
        found_existing=\$?

        # If the alias exists, but points to a different digest, then remove
        # the previous

        if [ "\$found_existing" -eq '0' ] && [ "\$existing_alias_digest" != "\$digest" ]; then
            echo "Removing alias \$alias on  \$existing_alias_digest to re-point it at \$digest"
            refgenie alias remove -a \$alias -d \$existing_alias_digest
        fi

        # The the alias didn't already exist, or did not match the current
        # digest (and so was removed), then add it

        if [ "\$found_existing" -eq '1' ] || [ "\$existing_alias_digest" != "\$digest" ]; then
            echo "Aliasing \$digest to \$alias"
            refgenie alias set --aliases \$alias --digest \$digest
            if [ \$? -ne 0 ]; then
                echo "Aliasing \$assembly to \$alias failed" 1>&2
                exit 1
            fi
        else
            echo "\$alias was already pointing at \$digest"
        fi
    done
    """
}

