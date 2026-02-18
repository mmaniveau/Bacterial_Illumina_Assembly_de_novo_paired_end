#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Paramètres configurables du pipeline 
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

/*
params.sample_csv        --> Fichier CSV : en première colonne le sample_id, suivie d'une colonne R1 et d'une colonne R2 avec un séparateur "," UFT8 sans BOM 
params.fastp_q           --> Seuil de qualité Phred pour Fastp
params.fastp_u           --> Pourcentage maximal de bases de mauvaise qualité par read pour Fastp
params.min_contig_length --> Taille minimale d’un contig à conserver après filtrage avec seqtk -L
params.pilon_nb_iter     --> Nombre d’itérations de correction du génome avec Pilon
params.spades_mode       --> Permet de choisir le mode d'assemblage : careful, meta ou default

Ces paramètres peuvent être modifiés par l’utilisateur pour ajuster le pipeline sans toucher au code tel que : 

--sample_csv <chemin vers le CSV>, ex: /home/User/Assembly/Test.csv
--fastp_q <int> ≥ 0 ; ≤ 40
--fastp_u <int> ≥ 0 ; ≤ 100
--min_contig_length <int> ≥ 100
--pilon_nb_iter <int> ≥ 1
--spades_mode <careful ou meta ou default>

*/

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Test de validation des paramètres
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


def validate_int_params(param_map) {
    def errors = []

    param_map.each { param_name, param_config ->
        def value = param_config.value
        def type = param_config.type
        def minValue = param_config.min
        def maxValue = param_config.max

        if (type == 'int' && value != null) {
            if (!(value instanceof Integer)) {
                errors << "'$param_name' doit être un entier."
                return
            }

            if (minValue != null && value < minValue) errors << "'$param_name' doit être ≥ ${minValue}."
            if (maxValue != null && value > maxValue) errors << "'$param_name' doit être ≤ ${maxValue}."
        }
    }

    if (errors) {
        error "[ERROR] Paramètres invalides :\n" + errors.join("\n")
    }
}

def validate_spades_mode(mode) {
    // 1-Liste des valeurs autorisées 
    def allowed_modes = ['default', 'careful', 'meta']
    
    if (!allowed_modes.contains(mode)) {
        error "[ERROR] Paramètre spades_mode invalide : '${mode}'. Choisir parmi : ${allowed_modes.join(', ')}."
    }

    // 2-Mapping vers l'option SPAdes avec tirets
    def spades_flag = ''
    if (mode == 'careful') spades_flag = '--careful'
    else if (mode == 'meta') spades_flag = '--meta'
    else if (mode == 'default') spades_flag = ''

    return spades_flag
}

validate_int_params([
    fastp_q:           [value: params.fastp_q, type: 'int', min: 0, max: 40],
    fastp_u:           [value: params.fastp_u, type: 'int', min: 0, max: 100],
    min_contig_length: [value: params.min_contig_length, type: 'int', min: 100],
    pilon_nb_iter:     [value: params.pilon_nb_iter, type: 'int', min: 1]
])

validate_spades_mode(params.spades_mode ?: 'normal')

println "\nPipeline d’assemblage bactérien – démarrage"
println "========================================"
println "fastp_q           : ${params.fastp_q}"
println "fastp_u           : ${params.fastp_u}"
println "min_contig_length : ${params.min_contig_length}"
println "pilon_nb_iter     : ${params.pilon_nb_iter}"
println "spades_mode       : ${params.spades_mode}"
println "========================================\n"

/*
Pipeline Illumina – Assemblage bactérien complet en paired-end
FastQC --> fastp --> FastQC --> SPAdes --> QUAST --> Filtrage contigs --> QUAST --> Pilon ×3 itérations --> QUAST --> Prokka
*/

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// 1-FastQC initial des données brutes
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

/* 
Entrée : 
Les channels reads contenant les fastq bruts sont transmis via le workflow dans le process fastqc_trimmed
ex : raw_reads = tuple("sample1", "sample1_R1.fastq", "sample1_R2.fastq")

Traitement :
FastQC effectue un contrôle qualité sur les reads brutes afin de vérifier la qualité du séquençage

Sortie :
Les fichier sont générés dans le dossier work et copiés dans :
qc_fastqc_raw_reads/${sample_id}/

Le process sera réutilisé ultérieurement sur les reads après trimming

*/

process fastqc_eval {

    tag "$sample_id"
    input:
    tuple val(sample_id), path(r1), path(r2)
    val tag_directory

    output:
    path "*_fastqc.{zip,html}"

    publishDir "qc_fastqc_${tag_directory}/${sample_id}", mode: 'copy'

    script:
    """
    fastqc $r1 $r2
    """
}


// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// 2-Trimming fastp
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//Alternative possible : Trimmomatic

/* 

Entrée : 
Les channels reads contenant les fichiers fastq bruts sont transmis via le workflow dans le process trim_reads
ex : reads = tuple("sample1", "sample1_R1.fastq", "sample1_R2.fastq")

Traitement :
fastp effectue le trimming des reads pour améliorer la qualité des séquences :
-q X : bases avec un score Phred inférieur à X sont éliminées
-u Y : reads avec plus de Y % de bases de qualité inférieure à X sont supprimés

Sortie :
Les fichier sont générés dans le dossier work et copiés dans :
qc_fastqc_trim_reads/${sample_id}/

*/

process trim_reads {

    tag "$sample_id"
    input:
    tuple val(sample_id), path(r1), path(r2)

    output:
    tuple val(sample_id),
          path("${sample_id}_R1_trim.fastq"),
          path("${sample_id}_R2_trim.fastq")

    publishDir "trim_reads/${sample_id}", mode: 'copy'

    script:
    """
    fastp \
      -i $r1 \
      -I $r2 \
      -o ${sample_id}_R1_trim.fastq \
      -O ${sample_id}_R2_trim.fastq \
      -h ${sample_id}_fastp.html \
      -j ${sample_id}_fastp.json \
      -q ${params.fastp_q} \
      -u ${params.fastp_u} 
    """
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// 3-FastQC après trimming
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

/* 

Entrée : 
Les channels trimmed contenant les fichiers fastq trimmés sont transmis via le workflow dans le process fastqc_trimmed
ex : trimmed = tuple("sample1", "sample1_R1_trim.fastq", "sample1_R2_trim.fastq")

Traitement :
FastQC effectue un contrôle qualité sur les reads après trimming, 
pour vérifier l’amélioration de la qualité des bases, l'élimination d'adapteurs etc.

Sortie :
Les fichier sont générés dans le dossier work et copiés dans :
qc_fastqc_trim_reads/${sample_id}/

*/


// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// 4-Assemblage SPAdes
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//Alternative possible : MEGAHIT, SKESA

/* 
Entrée : 
Les channels trimmed contenant les fastq trimmés vont être transmis via le workflow dans le process spades_assembly
ex : trimmed = tuple("sample1", "sample1_R1_trim.fastq", "sample1_R2_trim.fastq")

Traitement :
SPAdes assemble les reads trimmés en contigs pour chaque échantillon.
- SPAdes crée un sous-dossier par échantillon : spades_${sample_id}/
- Le fichier contigs.fasta est produit dans ce sous-dossier, ex :

  .../work/batch_number/hash_process/spades_sample1/contigs.fasta

Pour que Nextflow puisse détecter contigs.fasta comme output du process :

cp spades_${sample_id}/contigs.fasta ${sample_id}_contigs.fasta

copie le fichier "contigs.fasta" -hors du sous-dossier SPAdes-, pour le placer àa la racine dans : 

  .../work/batch_number/hash_process/contigs.fasta

Ceci est nécessaire car Nextflow ne détecte que les fichiers déclarés dans l'output au niveau du chemin : spades_sample1\contigs.fasta
Si le fichier reste dans le sous-dossier, il ne sera pas retrouvé pour la suite du pipeline.

Sortie :
Les fichier sont générés dans le dossier work et copiés dans :
spades/${sample_id}/${sample_id}_contigs.fasta

*/

spades_mode_with_flag = validate_spades_mode(params.spades_mode)

process spades_assembly {

    tag "$sample_id"
    input:
    tuple val(sample_id), path(r1), path(r2)

    output:
    tuple val(sample_id), path("${sample_id}_contigs.fasta")

    publishDir "spades/${sample_id}", mode: 'copy'

    script:
    """
    spades.py -1 $r1 -2 $r2 -o spades_${sample_id} $spades_mode_with_flag

    #Spades donne des noms de contig très long cette commande
    awk '/^>/ {print ">C" ++i; next} {print}' spades_${sample_id}/contigs.fasta > ${sample_id}_contigs.fasta
    """
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// 5-QUAST
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/* 
Entrée : 
Le process quast_eval reçoit un tuple contenant l'identifiant de l'échantillon et le chemin vers le fichier à analyser :
ex : contigs = tuple("sample1", "sample1_contigs.fasta")
Le process reçoit également un argument "tag_directory" (val) qui permet de personnaliser le nom du dossier de sortie
(ex. "spades", "contig_filtered", "pilon").

Traitement :
Cette étape permet d'évaluer la qualité de l'assemblage du génome (contigs) à différentes étapes du pipeline :
- après assemblage SPAdes
- après filtrage des contigs
- après correction Pilon
QUAST utilisera différents critères d'évaluation : nombre de contigs, N50, taille totale, etc.

Sortie :
Les fichier sont générés dans le dossier work et copiés dans :
  quast_${sample_id}_${tag_directory}

*/

process quast_eval {

    tag "$sample_id"
    input:
    tuple val(sample_id), path(final_genome)
    val tag_directory

    output:
    path "quast_${sample_id}_${tag_directory}"

    publishDir "quast/${sample_id}_${tag_directory}", mode: 'copy'

    script:
    """
    mkdir -p quast_${sample_id}_${tag_directory}
    quast.py $final_genome -o quast_${sample_id}_${tag_directory}
    """
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// 6-Filtrage contigs 
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

/* 
Entrée : 
Les channels contigs contenant les reads assemblés en contigs vont être transmis via le workflow dans le process filter_contigs 
ex : contigs = tuple("sample1", "sample1_contigs.fasta")

Traitement :
Cette étape permet de filtrer les contigs selon leur longueur minimale définie par l'option "-L X" de "seqtk seq"
- Les contigs courts sont souvent mal assemblés, contiennent de nombreuses séquences répétées et sont difficilement annotables
- Seuls les contigs d'une taille ≥ X seront conservés

Sortie :
Les fichiers filtrés sont générés dans le dossier work et copiés dans :
      contigs_filtered/${sample_id}/contigs_filtered.fasta
-->Ce fichier est déclaré dans l’output du process et utilisé pour les étapes suivantes du pipeline

*/

process filter_contigs {

    tag "$sample_id"
    input:
    tuple val(sample_id), path(raw_contigs)

    output:
    tuple val(sample_id), path("${sample_id}_contigs_filtered.fasta")

    publishDir "contigs_filtered/${sample_id}", mode: 'copy'

    script:
    """
    seqtk seq -L $params.min_contig_length $raw_contigs > ${sample_id}_contigs_filtered.fasta
    """
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// 7-QUAST
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// 8-Pilon
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

/* 
Entrée : 
Les channels contenant :
- les contigs filtrés : filtered_contigs = tuple("sample1", "sample1_contigs_filtered.fasta")
- les reads trimmés : reads_trimmed = tuple("sample1", "sample1_R1_trim.fastq", "sample1_R2_trim.fastq")
sont transmis via le workflow dans le process pilon_loop
sous la forme combinée (pilon_input) tel que:
  pilon_input = tuple("sample1", [r1, r2], contigs_filtered)
où r1_trim et r2_trim sont convertis en tableau : [r1, r2]


Traitement : 
À chaque itération, Pilon produit un fichier corrigé :
      .../work/batch_number/hash_process/pilon_iter_i/pilon.fasta
  qui contient le génome corrigé pour l'itération i.

Ce fichier est ensuite copié dans :
      .../work/batch_number/hash_process/pilon_work.fasta
  à la racine work du process.
  --> Il écrase la version précédente.
  --> Il sert de génome d'entrée pour l'itération suivante
  --> Il n'y a qu'un seul pilon_work.fasta à la fois dans ce dossier


Sortie :
Une copie nommée par échantillon et numéro d'itération est également conservée :
      pilon_iter_i/${sample_id}_pilon_iter_i.fasta

À la fin des X itérations, le fichier final est copié et renommé :
      ${sample_id}_pilon_final.fasta
  --> C'est ce fichier qui est déclaré dans l'output et récupéré par Nextflow pour les étapes suivantes.
*/

process pilon_loop {

    tag "$sample_id"

    input:
    tuple val(sample_id), path(reads_trimmed), path(filtered_contigs)

    output:
    tuple val(sample_id), path("pilon_iter_*/*.fasta"), path("${sample_id}_pilon_final.fasta")

    publishDir "pilon/${sample_id}", mode: 'copy'

/* Le fichier filtered_contigs est d'abord copié dans le répertoire work de pilon pour ne pas risquer de l'écraser ailleurs

bwa index : Indexation des contigs (génome de réference ici)

bwa mem : Alignement des reads R1/R2 trimmés ${reads[0]} ${reads[1]} contenus dans un tableau sur les contigs
à partir du tuple pilon_input résultant de la jointure des tuples contig et trimmed sur clé sample_id
    -->samtools view : conversion du flux l'alignement .sam en .bam
    -->samtools sort : Tri du flux d'alignement .bam

samtools index : Indexation du fichier d'alignement .bam en .bai
*/

    script:
    """
    #!/bin/bash
    set -euo pipefail

    cp $filtered_contigs pilon_work.fasta

    for i in \$(seq 1 $params.pilon_nb_iter); do
        # Nettoyage des index BWA et BAM de l'itération précédente
        rm -f align_sorted.bam align_sorted.bam.bai pilon_work.fasta.*

        bwa index pilon_work.fasta

        bwa mem pilon_work.fasta ${reads_trimmed[0]} ${reads_trimmed[1]} | \
        samtools view -b - | \
        samtools sort -o align_sorted.bam

        samtools index align_sorted.bam

        pilon \
          --genome pilon_work.fasta \
          --frags align_sorted.bam \
          --outdir pilon_iter_\${i} \
          --output ${sample_id}_pilon_iter_\${i}

        sed -i 's/_pilon//g' pilon_iter_\${i}/${sample_id}_pilon_iter_\${i}.fasta
        cp pilon_iter_\${i}/${sample_id}_pilon_iter_\${i}.fasta pilon_work.fasta
    done

    cp pilon_work.fasta ${sample_id}_pilon_final.fasta
    
    # Nettoyage final des index BWA
    rm -f pilon_work.fasta.*
    """
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// 9-QUAST
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// 10-Prokka
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

/* 
Entrée : 
Les channels filtered_contigs contenant les reads assemblés en contigs et polis après correction par Pilon
vont être transmis via le workflow dans le process prokka_annot
ex : pilon_final = tuple("sample1", "sample1_pilon_final.fasta")

Traitement : 
Commande awk : permet de renommer les contigs dans les fichiers de manière plus courte 
--> si les noms de contigs sont trop longs l'annotation Prokka échoue

Sortie :
Dossier : prokka/${sample_id}/ contenant tous les fichiers d’annotation

*/

process prokka_annot {

    tag "$sample_id"

    input:
    tuple val(sample_id), path(genome)

    output:
    path "prokka_${sample_id}"

    publishDir "prokka/${sample_id}", mode: 'copy'

    script:
    """
    # Renommer les contigs pour Prokka (IDs courts)
    awk '/^>/ {print ">contig" ++i; next} {print}' $genome > genome_renamed.fasta

    prokka \
      --outdir prokka_${sample_id} \
      --prefix ${sample_id} \
      --compliant \
      genome_renamed.fasta
    """
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// WORKFLOW
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

/* 
Le workflow fonctionne comme un //main() en Python avec appel de fonctions (process ici)
Il appelle des process en consommant des channels et en produisant de nouveaux channels
Les channels sont des flux de données permettant de transférer des informations d'un process à l'autre
Les tuples dans le script ne contiennent pas directement des fichiers mais des pointeurs permettant de remonter les chemins des fichiers
*/


workflow {

    Raw_Reads_QC()

    Trimming_Fastp(
        Raw_Reads_QC.out.raw_reads
    )

    Assembly_Spades(
        Trimming_Fastp.out.trimmed_reads
    )

    Filter_Contigs(
        Trimming_Fastp.out.trimmed_reads,
        Assembly_Spades.out.raw_contigs
    )

    Polishing_Pilon(
        Trimming_Fastp.out.trimmed_reads,
        Filter_Contigs.out.filtered_contigs
    )

    Annotation_Prokka(
        Polishing_Pilon.out.prokka_input
    )
}


workflow Raw_Reads_QC {

    main:

        if (!params.sample_csv) {
            error "[ERROR] Vous devez fournir un CSV avec --sample_csv <chemin_fichier.csv>"
        }

        def csvFile = file(params.sample_csv)
        if (!csvFile.exists()) {
            error "[ERROR] Fichier CSV '${params.sample_csv}' introuvable."
        }

        // Lire le CSV avec splitCsv 
        def samples = Channel
            .fromPath(csvFile)
            .splitCsv(header:true)

        // Construire un channel simple (id, r1, r2)
        raw_reads = samples.map { row ->
            def sample_id = row.sample_id.trim()
            def r1 = file(row.R1.trim())
            def r2 = file(row.R2.trim())
            
            // Validation de l'existence des fichiers pour éviter un crash en cours de route
            if (!r1.exists()) {
                error "[ERROR] Fichier R1 introuvable pour '${sample_id}': ${r1}"
            }
            if (!r2.exists()) {
                error "[ERROR] Fichier R2 introuvable pour '${sample_id}': ${r2}"
            }
            tuple(sample_id, r1, r2)
        }

        // Affiche les samples
        raw_reads.subscribe { id, r1, r2 ->
            println "[INFO] Sample : ${id}, R1=${r1.name}, R2=${r2.name}"
        }

        // Appel du process FastQC
        fastqc_eval(raw_reads, "raw")

    emit:
        raw_reads
}


workflow Trimming_Fastp {

    take:
        raw_reads

    main:
        trimmed_reads = trim_reads(raw_reads)

        fastqc_eval(trimmed_reads, "trimmed")

    emit:
        trimmed_reads
}

workflow Assembly_Spades {

    take:
        trimmed_reads

    main:
        raw_contigs = spades_assembly(trimmed_reads)

        quast_eval(raw_contigs, "spades")

    emit:
        raw_contigs
}

workflow Filter_Contigs {

    take:
        trimmed_reads
        raw_contigs

// Filtrage des contigs courts

    main:
        filtered_contigs = filter_contigs(raw_contigs)

//Analyse qualité avec Quast sur les contigs après filtrage pour vérifier si trop ou pas assez de contigs ont été éliminés

        quast_eval(filtered_contigs, "contig_filtered")

    emit: 
        filtered_contigs
}


workflow Polishing_Pilon {

    take:
        trimmed_reads
        filtered_contigs

/*
Jointure des tuples pour préparer l'input de Pilon :
filtered = tuple(sample_id, contigs) +
trimmed = tuple(sample_id, trim_r1, trim_r2)

Après la jointure "join" sur la clé "sample_id" :

pilon_input = tuple(sample_id, contigs_filtered, trim_r1, trim_r2)

r1_trim, r2_trim dans le tuple input sont ensuite convertis en un objet de type tableau de sorte que 
pilon_final = tuple(sample_id, [r1,r2], contigs) 
où [r1,r2] pourront être appelés dans pilon_loop par ${reads_trimmed[0]} et ${reads_trimmed[1]}

*/
    main:
        pilon_input = filtered_contigs.join(trimmed_reads)

        pilon_out = pilon_loop(
            pilon_input.map { id, contigs, r1, r2 ->
                tuple(id, [r1, r2], contigs)
            }
        )

        
//Analyse qualité avec Quast sur le génome final corrigé

        quast_eval(
            pilon_out.map { id, iters, final_fasta -> tuple(id, final_fasta) },
            "pilon"

        )
        
        prokka_input = pilon_out.map { id, iters, final_fasta ->
            tuple(id, final_fasta)        
            }
        
    emit: 
        prokka_input

}

workflow Annotation_Prokka {

    take:
        prokka_input

// Annotation finale
    main:
        prokka_annot(prokka_input)

}




