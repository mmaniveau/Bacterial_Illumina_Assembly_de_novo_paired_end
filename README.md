# Pipeline d’Assemblage Bactérien Illumina _de novo_ (Paired-End)

## Description
⚠️ Ce pipeline a été développé **à titre pédagogique**, avec de nombreux commentaires pour expliquer chaque étape.
Il est conçu pour l’assemblage _de novo_ de génomes bactériens à partir de reads paired-end Illumina, et est voué à évoluer pour intégrer d’autres commentaires supplémentaires et options, comme le traitement de reads single-end.


Ce pipeline Nextflow permet l’assemblage complet de génomes bactériens à partir de reads Illumina paired-end.  
Il inclut toutes les étapes nécessaires pour produire un génome poli et annoté :

1. Contrôle qualité initial des reads avec **FastQC**
2. Trimming des reads avec **fastp**
3. Contrôle qualité post-trimming avec FastQC
4. Assemblage du génome avec **SPAdes**
5. Évaluation de l’assemblage avec **QUAST**
6. Filtrage des contigs courts avec **Seqtk**
7. Polissage du génome avec **Pilon**, sur un nombre d’itérations configurable
8. Annotation finale du génome avec **Prokka**

Le pipeline est entièrement modulable et peut être appliqué sur différents datasets sans modification du code.

---

## Prérequis
Les outils suivants doivent être installés et accessibles dans le `$PATH` :

- Nextflow
- FastQC
- fastp
- SPAdes
- QUAST
- Seqtk
- Pilon
- BWA
- SAMtools
- Prokka
- awk (optionnel, utilisé pour renommer les contigs avant Prokka)

---

## Installation et emplacement des reads

Placez vos fichiers **FASTQ paired-end** dans un dossier.  

- Si aucun dossier n’est indiqué via `--reads_dir`, le pipeline utilise le **dossier courant (`./`)** et traite tous les fichiers FASTQ qu’il y trouve.  
- Cela permet de lancer le pipeline directement sur les fichiers présents dans le répertoire où Nextflow est exécuté.

### Détection automatique des paires R1/R2

Le workflow **Raw_Reads_QC** détecte automatiquement toutes les paires de reads FASTQ dans le dossier fourni.  

**Fonctionnement :**

1. Chaque fichier FASTQ est analysé pour identifier un **suffixe indiquant R1 ou R2**.
2. Suffixes reconnus pour R1 : `R1`, `1`, `forward`, `fw`, `F`, `f`  
   Suffixes reconnus pour R2 : `R2`, `2`, `reverse`, `rev`, `R`, `r`
3. Tout texte après le suffixe (par exemple `_001.fastq`) est considéré comme un suffixe additionnel et **remis derrière le nom de base** de l’échantillon pour générer un identifiant unique.

**Exemples :**

- `SampleA_R1.fastq` et `SampleA_R2.fastq` → `sample_id = SampleA`
- `SampleB-fw_001.fastq` et `SampleB-rev_001.fastq` → `sample_id = SampleB_001`
- `my-sample.F` et `my-sample.R` → `sample_id = my-sample`

Le pipeline crée ensuite un **tuple** pour chaque échantillon contenant : `(sample_id, path_R1, path_R2)`  

Si une paire est incomplète, le pipeline s'arrête avec une erreur.

### Sortie affichée pour chaque sample détecté :
```

[INFO] Sample détecté : SampleA
R1: /chemin/vers/SampleA_R1.fastq
R2: /chemin/vers/SampleA_R2.fastq

````

Ainsi, quelle que soit la variante de suffixe utilisée pour indiquer R1/R2, le pipeline **reconstruit correctement l'identifiant de l'échantillon et lie les fichiers correspondants**.

---

## Paramètres configurables

Le pipeline peut être configuré via la ligne de commande Nextflow. Les paramètres modifiables sont :

- `fastp_q` : seuil de qualité Phred pour fastp (défaut 20, 0–40)
- `fastp_u` : pourcentage maximal de bases de mauvaise qualité par read (défaut 20, 0–100)
- `min_contig_length` : longueur minimale d’un contig à conserver (défaut 500, ≥100)
- `pilon_nb_iter` : nombre d’itérations pour Pilon (défaut 3, ≥1)
- `reads_dir` : dossier contenant les reads (défaut `./`)
- `spades_mode` : mode d’assemblage SPAdes (`default`, `careful`, `meta`; défaut `default`)
- `reads_dir` : chemin vers les fichiers à analyser

**Exemple de lancement avec paramètres personnalisés :**
```bash
nextflow run main.nf --reads_dir ./my_reads --fastp_q 25 --fastp_u 15 --min_contig_length 1000 --pilon_nb_iter 2 --spades_mode careful 
````

---

## Sorties

Pour chaque échantillon, le pipeline produit :

* **FastQC** : rapports HTML et ZIP pour reads bruts et trimmés
* **SPAdes** : contigs assemblés (`*_contigs.fasta`)
* **QUAST** : statistiques d’assemblage à différentes étapes
* **Seqtk** : contigs filtrés (`*_contigs_filtered.fasta`)
* **Pilon** : génomes polis après chaque itération (`*_pilon_iter_*.fasta`) et final (`*_pilon_final.fasta`)
* **Prokka** : annotation complète dans `prokka/${sample_id}`

---

## Organisation des dossiers de sortie

* FastQC initial : `qc_fastqc_raw/${sample_id}/`
* Trimming : `trim_reads/${sample_id}/`
* FastQC post-trimming : `qc_fastqc_trimmed/${sample_id}/`
* Assemblage SPAdes : `spades/${sample_id}/`
* Contigs filtrés : `contigs_filtered/${sample_id}/`
* Polissage Pilon : `pilon/${sample_id}/`
* Annotation Prokka : `prokka/${sample_id}/`

---

## Exécution du pipeline

Pour lancer le pipeline :

```bash
nextflow run main.nf
```

## Validation des entrées

Le pipeline effectue plusieurs vérifications avant de commencer :

1. Le dossier de reads doit exister.
2. Chaque R1 doit avoir un R2 correspondant.
3. Les paramètres numériques (`fastp_q`, `fastp_u`, `min_contig_length`, `pilon_nb_iter`) doivent être valides.
4. Le mode SPAdes doit être l’un des suivants : `default`, `careful`, `meta`.

---

## Schéma du pipeline

```
Raw Reads (FASTQ)
      │
      ▼
   FastQC
      │
      ▼
   fastp Trimming
      │
      ▼
   FastQC (trimmed)
      │
      ▼
    SPAdes Assembly
      │
      ▼
     QUAST
      │
      ▼
  Filter Contigs (seqtk)
      │
      ▼
     QUAST
      │
      ▼
   Pilon Polishing (xN iter)
      │
      ▼
     QUAST
      │
      ▼
    Prokka Annotation
      │
      ▼
  Annotated Genome (GFF, GBK, etc.)
```

---

## Remarques

* Conçu uniquement pour **reads paired-end Illumina**.
* Les étapes peuvent être remplacées par d’autres outils (ex. Trimmomatic pour fastp, MEGAHIT pour SPAdes).
* Prokka renomme automatiquement les contigs si nécessaire pour éviter des erreurs liées à des noms trop longs.
* Si aucun dossier n’est indiqué via `--reads_dir`, le pipeline s’exécute sur tous les fichiers FASTQ présents dans le répertoire courant.

---
