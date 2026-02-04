# Pipeline d’Assemblage Bactérien Illumina _de novo_ (Paired-End)

## Description
⚠️ Ce pipeline a été développé **à titre pédagogique**, avec de nombreux commentaires pour expliquer chaque étape.
Il est conçu pour l’assemblage _de novo_ de génomes bactériens à partir de reads paired-end Illumina, et est voué à évoluer pour intégrer d’autres commentaires supplémentaires et options, comme le traitement de reads single-end ou en collection.


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


## Installation et fourniture des reads

Le pipeline utilise désormais un **fichier CSV** pour définir explicitement les reads **FASTQ paired-end** à analyser.

### Fichier CSV des échantillons

Les fichiers FASTQ doivent être déclarés dans un **CSV fourni via l’option `--sample_csv`**.

Le CSV doit contenir au minimum les colonnes suivantes (les noms sont **case-sensitifs**) :


sample_id,R1,R2

Des colonnes supplémentaires (ex. Meta) peuvent être présentes mais ne sont pas utilisées par le pipeline principal.

Exemple de CSV valide
sample_id,R1,R2,meta
Sample_1,/home/user/data/Sample_1_R1.fastq,/home/user/data/Sample_1_R2.fastq,Sel


sample_id : identifiant unique de l’échantillon

R1 : chemin complet vers le fichier FASTQ R1

R2 : chemin complet vers le fichier FASTQ R2


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
- `sample_csv` : chemin vers le csv listant les fichiers à analyser avec leur chemin 

**Exemple de lancement avec paramètres personnalisés :**
```bash
nextflow run main.nf --sample_csv ./my_reads.csv --fastp_q 25 --fastp_u 15 --min_contig_length 1000 --pilon_nb_iter 2 --spades_mode careful 
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
---
