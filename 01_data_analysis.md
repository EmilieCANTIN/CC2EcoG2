Analyse\_des\_données\_de\_la\_rade\_de\_Brest
================

  - [Méthodes](#méthodes)
      - [Des lectures aux tableaux](#des-lectures-aux-tableaux)
      - [Filtrer les données](#filtrer-les-données)
      - [Filtrez les lectures en amont et en
        aval](#filtrez-les-lectures-en-amont-et-en-aval)
  - [Connaître les taux d’erreur](#connaître-les-taux-derreur)
      - [Exemple d’inférence](#exemple-dinférence)
  - [Fusions des lectures et élimination des
    chimères](#fusions-des-lectures-et-élimination-des-chimères)
      - [Construire un tableau
        séquentiel](#construire-un-tableau-séquentiel)
      - [Supprimer les chimères](#supprimer-les-chimères)
  - [“Track reads through the
    pipeline”](#track-reads-through-the-pipeline)
  - [Attribuer une taxonomie](#attribuer-une-taxonomie)

# Méthodes

## Des lectures aux tableaux

1/ quelles sont les influences relatives de la profondeur et de la
saison sur la structure des communautes planctoniques de la rade de
Brest ? 2/ Quels sont les biomarqueurs de saison (hiver et ete) ?

``` r
library("dada2")
```

    ## Loading required package: Rcpp

    ## Warning: multiple methods tables found for 'which'

Ce premier code permet d’importer les données de l’étude de la rade de
Brest, à partir d’un ensemble de fichiers fastq.Ici, on définit une
variable chemin path, afin de pouvoir accéder à ces données.

``` r
path <- "~/CC2EcoG2/Seqreunies" # MODIFIER le répertoire contenant les fichiers fastq après la décompression
list.files(path)
```

    ##  [1] "filtered"                            "Station5_Fond1_10sept14_R1.fastq"   
    ##  [3] "Station5_Fond1_10sept14_R2.fastq"    "Station5_Fond1_11mars15_R1.fastq"   
    ##  [5] "Station5_Fond1_11mars15_R2.fastq"    "Station5_Fond2_10sept14_R1.fastq"   
    ##  [7] "Station5_Fond2_10sept14_R2.fastq"    "Station5_Fond2_11mars15_R1.fastq"   
    ##  [9] "Station5_Fond2_11mars15_R2.fastq"    "Station5_Fond3_10sept14_R1.fastq"   
    ## [11] "Station5_Fond3_10sept14_R2.fastq"    "Station5_Median1_10sept14_R1.fastq" 
    ## [13] "Station5_Median1_10sept14_R2.fastq"  "Station5_Median2_10sept14_R1.fastq" 
    ## [15] "Station5_Median2_10sept14_R2.fastq"  "Station5_Surface1_10sept14_R1.fastq"
    ## [17] "Station5_Surface1_10sept14_R2.fastq" "Station5_Surface1_11mars15_R1.fastq"
    ## [19] "Station5_Surface1_11mars15_R2.fastq" "Station5_Surface2_10sept14_R1.fastq"
    ## [21] "Station5_Surface2_10sept14_R2.fastq" "Station5_Surface2_11mars15_R1.fastq"
    ## [23] "Station5_Surface2_11mars15_R2.fastq"

## Filtrer les données

On filtre les séquences de faible qualité, puis on les enlève. On
demande ici d’afficher les “moins bons”.

``` r
# Le tri permet de s'assurer que les lectures en avant et en arrière sont dans le même ordre
fnFs <- sort(list.files(path, pattern="_R1.fastq"))
fnRs <- sort(list.files(path, pattern="_R2.fastq"))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "R"), `[`, 1)
# Specify the full path to the fnFs and fnRs
fnFs <- file.path(path, fnFs)
fnRs <- file.path(path, fnRs)
```

On sait que plus on se rapproche de la fin des séquençages, moins bonne
sera leur qualité. En effet, on remarque que pour les lectures avant
(deux premiers graphes), le score de qualité moyen ne descend jamais en
dessous de 30. Au contraire, les graphes incarnant la fin des lectures
montrent un score de qualité plus bas (\~25). Ce type de chiffre
représente la probabilité que ce ne soit pas le bon nucléotide
d’appelé. De ce fait, avec un Q30 en début de séquences, il y a une
chance sur 1000 que ce soit le cas.

``` r
library(dada2)
plotQualityProfile(fnFs[1:2])
```

![](01_data_analysis_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

``` r
plotQualityProfile(fnRs[1:2])
```

![](01_data_analysis_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

On voit bien ici la moins bonne qualité des fins de séquences, que ce
soit en été ou en hiver. En effet, les scores de qualités baissent vers
la position 240 pour les premières lectures, et plutôt vers la position
200 pour les lectures arrières. En prenant ces informations en compte,
on va pouvoir dans un premier temps créer des variables pour les
fichiers filtrés, puis appliquer la fonction filterAndTrim.

``` r
filt_path <- file.path(path, "filtered") # Placez les fichiers filtrés dans le sous-répertoire "filtered"
if(!file_test("-d", filt_path)) dir.create(filt_path)
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names
sample.names
```

    ##  [1] "Station5_Fond1_10sept14_"    "Station5_Fond1_11mars15_"   
    ##  [3] "Station5_Fond2_10sept14_"    "Station5_Fond2_11mars15_"   
    ##  [5] "Station5_Fond3_10sept14_"    "Station5_Median1_10sept14_" 
    ##  [7] "Station5_Median2_10sept14_"  "Station5_Surface1_10sept14_"
    ##  [9] "Station5_Surface1_11mars15_" "Station5_Surface2_10sept14_"
    ## [11] "Station5_Surface2_11mars15_"

## Filtrez les lectures en amont et en aval

Cette fonction se base sur des fichiers contenant les lectures coupées
ayant passées les filtres. “TrimLeft” permet de retirer les primers afin
de ne pas les intégrer aux séquences étudiées.

``` r
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,200), trimLeft=c(21,21),
              maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
              compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
head(out)
```

    ##                                    reads.in reads.out
    ## Station5_Fond1_10sept14_R1.fastq     159971    145448
    ## Station5_Fond1_11mars15_R1.fastq     175993    160423
    ## Station5_Fond2_10sept14_R1.fastq     197039    177018
    ## Station5_Fond2_11mars15_R1.fastq      87585     79989
    ## Station5_Fond3_10sept14_R1.fastq     117140    106150
    ## Station5_Median1_10sept14_R1.fastq   116519    106745

# Connaître les taux d’erreur

Ci-dessous, la fonction learnErrors permet d’estimer les taux d’erreurs
à partir d’un grand ensemble de données. Ainsi, les résultats ci-après
expriment le nombre de bases qui sera finalement utilisé, par rapport au
premier ensemble.

``` r
errF <- learnErrors(filtFs, multithread=TRUE)
```

    ## 105752691 total bases in 482889 reads from 3 samples will be used for learning the error rates.

``` r
errR <- learnErrors(filtRs, multithread=TRUE)
```

    ## 100755162 total bases in 562878 reads from 4 samples will be used for learning the error rates.

``` r
plotErrors(errF, nominalQ=TRUE)
```

    ## Warning: Transformation introduced infinite values in continuous y-axis
    
    ## Warning: Transformation introduced infinite values in continuous y-axis

![](01_data_analysis_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

``` r
library(ggplot2)
library(dada2)
plotErrors(errF, nominalQ=TRUE)
```

    ## Warning: Transformation introduced infinite values in continuous y-axis
    
    ## Warning: Transformation introduced infinite values in continuous y-axis

![](01_data_analysis_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

Les figures ci-dessus représentent les estimations des taux d’erreurs en
fonction du score de qualité. La ligne rouge incarne la tendance
générale du graphique. Ensuite, les points noirs reflètent le taux
d’erreurs observées, et la ligne noire le taux d’erreurs ajustées. On
peut donc observer ci-dessus la fréquence du taux d’erreur en fonction
du score de qualité. Aucune différence significative ne peut être
relevée entre errR et errF. En effet, on observe la même tendance :
moins il y a d’erreurs, plus le score de qualité augmente, ce qui est en
accord avec les résultats attendus.

## Exemple d’inférence

La fonction dada retire les erreurs de séquençage et renvoie la
composition déduite des échantillons.

``` r
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
```

    ## Sample 1 - 145448 reads in 37907 unique sequences.
    ## Sample 2 - 160423 reads in 35863 unique sequences.
    ## Sample 3 - 177018 reads in 47212 unique sequences.
    ## Sample 4 - 79989 reads in 20356 unique sequences.
    ## Sample 5 - 106150 reads in 30255 unique sequences.
    ## Sample 6 - 106745 reads in 28836 unique sequences.
    ## Sample 7 - 98823 reads in 25824 unique sequences.
    ## Sample 8 - 107427 reads in 26733 unique sequences.
    ## Sample 9 - 71082 reads in 17976 unique sequences.
    ## Sample 10 - 78645 reads in 20422 unique sequences.
    ## Sample 11 - 91534 reads in 24487 unique sequences.

``` r
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
```

    ## Sample 1 - 145448 reads in 45486 unique sequences.
    ## Sample 2 - 160423 reads in 41638 unique sequences.
    ## Sample 3 - 177018 reads in 55554 unique sequences.
    ## Sample 4 - 79989 reads in 23239 unique sequences.
    ## Sample 5 - 106150 reads in 34625 unique sequences.
    ## Sample 6 - 106745 reads in 31673 unique sequences.
    ## Sample 7 - 98823 reads in 29093 unique sequences.
    ## Sample 8 - 107427 reads in 28947 unique sequences.
    ## Sample 9 - 71082 reads in 21426 unique sequences.
    ## Sample 10 - 78645 reads in 22051 unique sequences.
    ## Sample 11 - 91534 reads in 28266 unique sequences.

# Fusions des lectures et élimination des chimères

``` r
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
```

    ## 117318 paired-reads (in 5196 unique pairings) successfully merged out of 141000 (in 21451 pairings) input.

    ## 138940 paired-reads (in 4296 unique pairings) successfully merged out of 156462 (in 15709 pairings) input.

    ## 142188 paired-reads (in 6989 unique pairings) successfully merged out of 171439 (in 27056 pairings) input.

    ## 67622 paired-reads (in 2721 unique pairings) successfully merged out of 77764 (in 9556 pairings) input.

    ## 83613 paired-reads (in 3458 unique pairings) successfully merged out of 102224 (in 16304 pairings) input.

    ## 86212 paired-reads (in 3348 unique pairings) successfully merged out of 103447 (in 14293 pairings) input.

    ## 80661 paired-reads (in 2727 unique pairings) successfully merged out of 95866 (in 12350 pairings) input.

    ## 89385 paired-reads (in 3073 unique pairings) successfully merged out of 104354 (in 12135 pairings) input.

    ## 59716 paired-reads (in 1939 unique pairings) successfully merged out of 68711 (in 7974 pairings) input.

    ## 66157 paired-reads (in 1763 unique pairings) successfully merged out of 76701 (in 8283 pairings) input.

    ## 75048 paired-reads (in 3149 unique pairings) successfully merged out of 88514 (in 12054 pairings) input.

``` r
# Inspecter le fichier "merger data.frame" du premier échantillon
head(mergers[[1]])
```

    ##                                                                                                                                                                                                                                                                                                                                                                                sequence
    ## 1     TACGAAGGGACCTAGCGTAGTTCGGAATTACTGGGCTTAAAGAGTTCGTAGGTGGTTGAAAAAGTTAGTGGTGAAATCCCAGAGCTTAACTCTGGAACTGCCATTAAAACTTTTCAGCTAGAGTATGATAGAGGAAAGCAGAATTTCTAGTGTAGAGGTGAAATTCGTAGATATTAGAAAGAATACCAATTGCGAAGGCAGCTTTCTGGATCATTACTGACACTGAGGAACGAAAGCATGGGTAGCGAAGAGGATTAGATACCCTCGTAGTCCATGCCGTAAACGATGTGTGTTAGACGTTGGAAATTTATTTTCAGTGTCGCAGGGAAACCGATAAACACACCGCCTGGGGAGTACGACCGCAAGGTT
    ## 2     TACGAAGGGACCTAGCGTAGTTCGGAATTACTGGGCTTAAAGAGTTCGTAGGTGGTTGAAAAAGTTGGTGGTGAAATCCCAGAGCTTAACTCTGGAACTGCCATCAAAACTTTTCAGCTAGAGTATGATAGAGGAAAGCAGAATTTCTAGTGTAGAGGTGAAATTCGTAGATATTAGAAAGAATACCAATTGCGAAGGCAGCTTTCTGGATCATTACTGACACTGAGGAACGAAAGCATGGGTAGCGAAGAGGATTAGATACCCTCGTAGTCCATGCCGTAAACGATGTGTGTTAGACGTTGGAAATTTATTTTCAGTGTCGCAGCGAAAGCGATAAACACACCGCCTGGGGAGTACGACCGCAAGGTT
    ## 3     TACGAAGGGACCTAGCGTAGTTCGGAATTACTGGGCTTAAAGAGTTCGTAGGTGGTTGAAAAAGTTGGTGGTGAAATCCCAGAGCTTAACTCTGGAACTGCCATCAAAACTTTTCAGCTAGAGTTTGATAGAGGAAAGCAGAATTTCTAGTGTAGAGGTGAAATTCGTAGATATTAGAAAGAATACCAATTGCGAAGGCAGCTTTCTGGATCATTACTGACACTGAGGAACGAAAGCATGGGTAGCGAAGAGGATTAGATACCCTCGTAGTCCATGCCGTAAACGATGTGTGTTAGACGTTGGAAATTTATTTTCAGTGTCGCAGCGAAAGCGATAAACACACCGCCTGGGGAGTACGACCGCAAGGTT
    ## 4     TACGAAGGGACCTAGCGTAGTTCGGAATTACTGGGCTTAAAGAGTTCGTAGGTGGTTGAAAAAGTTAGTGGTGAAATCCCAGAGCTTAACTCTGGAACTGCCATTAAAACTTTTCAGCTAGAGTATGATAGAGGAAAGCAGAATTTCTAGTGTAGAGGTGAAATTCGTAGATATTAGAAAGAATACCAATTGCGAAGGCAGCTTTCTGGATCATTACTGACACTGAGGAACGAAAGCATGGGTAGCGAAGAGGATTAGATACCCTCGTAGTCCATGCCGTAAACGATGTGTGTTAGACGTTGGAAATTTATTTTCAGTGTCGCAGCGAAAGCGATAAACACACCGCCTGGGGAGTACGACCGCAAGGTT
    ## 5     TACGAAGGGACCTAGCGTAGTTCGGAATTACTGGGCTTAAAGAGTTCGTAGGTGGTTGAAAAAGTTGGTGGTGAAATCCCAGAGCTTAACTCTGGAACTGCCATCAAAACTTTTCAGCTAGAGTATGATAGAGGAAAGCAGAATTTCTAGTGTAGAGGTGAAATTCGTAGATATTAGAAAGAATACCAATTGCGAAGGCAGCTTTCTGGATCATTACTGACACTGAGGAACGAAAGCATGGGTAGCGAAGAGGATTAGATACCCTCGTAGTCCATGCCGTAAACGATGTGTGTTAGACGTTGGAAATTTATTTTCAGTGTCGCAGGGAAACCGATAAACACACCGCCTGGGGAGTACGACCGCAAGGTT
    ## 6 TACGAGGGGTCCTAGCGTTGTCCGGATTTACTGGGCGTAAAGGGTACGTAGGCGTTTTAATAAGTTGTATGTTAAATATCTTAGCTTAACTAAGAAAGTGCATACAAAACTGTTAAGATAGAGTTTGAGAGAGGAACGCAGAATTCATGGTGGAGCGGTGACATGCGTAGATATCATGAGGAAAGTCAAATGCGAAGGCAGCCTTCTGGCTCAAAACTGACGCTGAGGTACGAAAGCGTGGGGAGCGAACAGGATTAGATACCCTGGTAGTCCACGCCGTAAACGATGAGTATTTGGTGCTGGGGGATTCGACCCTTTCAGTGCCGTAGCTAACGCGATAAATACTCCGCCTGGGGACTACGATCGCAAGATT
    ##   abundance forward reverse nmatch nmismatch nindel prefer accept
    ## 1      5170       1       2     29         0      0      2   TRUE
    ## 2      4129       2       1     29         0      0      2   TRUE
    ## 3      3757       3       1     29         0      0      2   TRUE
    ## 4      2481       1       1     29         0      0      2   TRUE
    ## 5      2182       2       2     29         0      0      2   TRUE
    ## 6      2132       5       9     25         0      0      1   TRUE

## Construire un tableau séquentiel

``` r
seqtab <- makeSequenceTable(mergers)
# Récupérer ou définir la dimension d'un objet
dim(seqtab)
```

    ## [1]    11 19426

``` r
# Contrôler la distribution des longueurs de séquence
table(nchar(getSequences(seqtab)))
```

    ## 
    ##  352  353  362  363  364  365  366  367  368  369  370  371  372  373  374  375 
    ##    1    1    1    1    4  183   27  165  184 5608 3594 2312 2613 2738  126 1770 
    ##  376  377  378  382  386 
    ##   90    4    1    1    2

Ici, il supprime les séquences reproduites en comparant chaque séquence
aux autres.

## Supprimer les chimères

``` r
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
```

    ## Identified 17869 bimeras out of 19426 input sequences.

``` r
dim(seqtab.nochim)
```

    ## [1]   11 1557

On peut donc dire que les chimères représentent environ 11% des
variantes de séquences.

``` r
sum(seqtab.nochim)/sum(seqtab)
```

    ## [1] 0.7769154

Il y a environ 23 pourcents de chimères.

# “Track reads through the pipeline”

``` r
library(dada2)
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# Si vous traitez un seul échantillon, supprimez les sapply : remplacez sapply(dadaFs, getN) par getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
```

    ##                             input filtered denoisedF denoisedR merged nonchim
    ## Station5_Fond1_10sept14_   159971   145448    142931    143292 117318   87962
    ## Station5_Fond1_11mars15_   175993   160423    158128    158473 138940  111552
    ## Station5_Fond2_10sept14_   197039   177018    173601    174591 142188  103668
    ## Station5_Fond2_11mars15_    87585    79989     78618     78926  67622   54711
    ## Station5_Fond3_10sept14_   117140   106150    103806    104338  83613   64259
    ## Station5_Median1_10sept14_ 116519   106745    104811    105173  86212   65559

# Attribuer une taxonomie

On va assigner une taxonomie aux données de la rade de Brest à partir de
ce qu’on pourra observer de Silva. Cela est rendu possible car l’ARN 16S
est un marqueur extrêmement bien classé. Tout d’abord, il nous faut
importer les données Silva, qui nous serviront à réaliser l’arbre
taxonomique. Vous pourrez retrouver les codes nécessaires dans
00\_data\_import.

``` r
taxa <- assignTaxonomy(seqtab.nochim, "~/silva_nr99_v138_train_set.fa.gz", multithread=TRUE)
```

``` r
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)
```

    ##      Kingdom    Phylum             Class                 Order            
    ## [1,] "Bacteria" "Proteobacteria"   "Alphaproteobacteria" "SAR11 clade"    
    ## [2,] "Bacteria" "Cyanobacteria"    "Cyanobacteriia"      "Synechococcales"
    ## [3,] "Bacteria" "Proteobacteria"   "Alphaproteobacteria" "SAR11 clade"    
    ## [4,] "Bacteria" "Proteobacteria"   "Alphaproteobacteria" "SAR11 clade"    
    ## [5,] "Bacteria" "Proteobacteria"   "Alphaproteobacteria" "SAR11 clade"    
    ## [6,] "Bacteria" "Actinobacteriota" "Acidimicrobiia"      "Actinomarinales"
    ##      Family             Genus                    
    ## [1,] "Clade I"          "Clade Ia"               
    ## [2,] "Cyanobiaceae"     "Synechococcus CC9902"   
    ## [3,] "Clade I"          "Clade Ia"               
    ## [4,] "Clade I"          "Clade Ia"               
    ## [5,] "Clade II"         NA                       
    ## [6,] "Actinomarinaceae" "Candidatus Actinomarina"

On peut ici remarquer la présence de plusieurs genres différents. On
pourra étudier leur répartition en fonction de la profondeur et de la
saison dans l’étude phyloseq.

``` r
save.image(file = "01_data_analysis_FinalEnv")
```
