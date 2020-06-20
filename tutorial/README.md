
![alt text](../markdown/mocca.png "")

# MOCCA
Copyright Bjørn Bredesen, 2011-2020

-------------------------------------------------

## Tutorial

In this tutorial, we demonstrate how to train models implemented in the Motif Occurrence Combinatorics Classification Algorithms (MOCCA) suite.
We previously presented the Support Vector Machine Motif Occurrence Combinatorics Classification Algorithm (SVM-MOCCA) (Bredesen *et al.* 2019), which we found to improve generalization to Polycomb/Trithorax Response Elements (PREs), a class of *cis*-regulatory elements (CREs) that maintains epigenetic memory.
SVM-MOCCA is a hierarchical method based on Support Vector Machines (SVMs) and motifs, where one SVM is trained per motif to classify its occurrences, with the feature space consisting of local dinucleotide and motif occurrence frequencies.
A key novelty of the MOCCA suite is the inclusion of the first polished and configurable implementation of SVM-MOCCA. Additionally, MOCCA includes a derivative method based on Random Forests (RFs), called the Random Forest Motif Occurrence Combinatorics Classification Algorithm (RF-MOCCA). Given their importance, we will focus on the training of these models first.
In addition, MOCCA supports the training of log-odds and classic SVM and RF models based on a variety of feature space formulations, which is treated in the final sections of the tutorial.


### Preparations
--------------------------

This tutorial assumes that MOCCA has been successfully installed. For directions regarding the installation of MOCCA, please see the main README file that accompanies MOCCA.

Sequence models in MOCCA make use of sequence motifs in order to learn distinguishing features of the different sequence classes. Accordingly, the motifs to be used must be specified. The simplest way of supplying MOCCA with a set of known motifs is by specifying them in IUPAC nucleotide codes (https://www.bioinformatics.org/sms/iupac.html) in an XML-file. We have supplied a sample XML-file, with the motifs that we used previously for genome-wide prediction of PREs (Bredesen *et al* 2019). The motifs are given in the file `motifs2019.xml`. Alternatively, motifs can be specified individually with commandline arguments (see `mocca --help`). In addition, motifs can be specified using Position Weight Matrices (PWMs). In order to demonstrate the use of PWMs, we additionally supply the XML file `motifs2019expho.xml`, which has the motifs for the DNA-binding factor *pho* removed, and we include a PWM-motif for *pho* (in the file `UmassPGFE_PSSM_pho_SOLEXA_5_20200302.pssm`) acquired from the Fly Factor Survey.

MOCCA uses binary or multi-class machine learning, and thus requires training data of at least two different classes. Classes and training data can be specified manually using commandline arguments. MOCCA also supplies a mode that automatically constructs negative training data based on a set of positives and an input genome, given in FASTA format &ndash; here referred to as "assisted training".

In this tutorial, we will use the *D. melanogaster* genome assembly R5, which can be downloaded from FlyBase:
 * **Either**: Download and unpack the genome sequence to the folder with the tutorial data: ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r5.57_FB2014_03/fasta/dmel-all-chromosome-r5.57.fasta.gz
 * **Or**: In a Linux terminal, navigate to the folder with the tutorial data, and execute the following command: `wget ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r5.57_FB2014_03/fasta/dmel-all-chromosome-r5.57.fasta.gz && gunzip dmel-all-chromosome-r5.57.fasta.gz`


### SVM-MOCCA &ndash; Assisted training
--------------------------

In this section, we will train SVM-MOCCA to predict Polycomb/Trithorax Response Elements (PREs) using assisted training and genome-wide Polycomb target sites from the Kahn *et al.* (2014) experimental set. With assisted learning, MOCCA requires only a positive training set and a genome, and MOCCA generates negative data using generative models. We extracted the Polycomb target site regions determined by Kahn *et al.* (2014) from their Supplementary Table 3, we resized all regions to a length of 3kb each, and we extracted the underlying sequences from the *D. melanogaster* genome release 5, using Gnocis. The resulting sequences are supplied in the file `KahnPcG.fa`.

With this set up, we are ready to train SVM-MOCCA. We will train SVM-MOCCA using the IUPAC motifs with *pho* excluded, and with a PWM motif for *pho*. In a terminal, navigate to the tutorial folder, and run:
```
mocca -auto:FASTA KahnPcG.fa -genome:FASTA dmel-all-chromosome-r5.57.fasta -motif:XML motifs2019expho.xml -motif:PWM UmassPGFE_PSSM_pho_SOLEXA_5_20200302.pssm 0 
```

MOCCA will then execute the following steps:
 * Load the input files.
 * Divide the training sequences (`KahnPcG.fa`) into independent training and test sets.
 * Train a 4th order Markov chain on the training sequences and generate negative training and test data.
    - This class of negatives is called "dummy CREs", due to being generated by a model trained on PREs.
    - The negative test data is 100 times larger than the positive test data, in order to reflect the expected imbalance in the genome.
    - Note: The order of the Markov chain can be changed, or an i.i.d. model can optionally be used (see `mocca --help`).
 * Train a 4th order Markov chain on the genome and generate negative training, test and calibration data.
    - This class of negatives is called "dummy genomic", due to being generated by a model trained genome-wide.
    - The calibration data is the size of the genome.
    - Note: The order of the Markov chain can be changed, or an i.i.d. model can optionally be used (see `mocca --help`).
 * Calibrate the PWM-threshold for one occurrence expected per kilobase.
 * Train an SVM-MOCCA model, with three classes (CRE, dummy CRE and dummy genomic).
 * Calibrate the threshold for an expected genome-wide precision of 80% (can be changed).
 * Print out validation statistics.
 * *Note*: This may take some time to execute.

We usually also want to save the predictions. In order to do so, at least one additional argument must be specified. In order to output genome-wide predictions as non-overlapping regions in the General Feature Format (https://www.ensembl.org/info/website/upload/gff.html), add `-predict:GFF predictions.gff`. To output all genome-wide window scores, add `-predict:Wig predictions.wig`.

```
mocca -auto:FASTA KahnPcG.fa -genome:FASTA dmel-all-chromosome-r5.57.fasta -motif:XML motifs2019expho.xml -motif:PWM UmassPGFE_PSSM_pho_SOLEXA_5_20200302.pssm 0 -predict:GFF predictions.gff -predict:Wig predictions.wig
```

After running this, MOCCA should have generated the files `predictions.gff` and `predictions.wig` in the same folder.

Note: If prediction is to be performed on only a subset of the chromosomes specified in the input genome FASTA file, the FASTA file must be filtered in advance of prediction.


### SVM-MOCCA &ndash; Core-CRE prediction
--------------------------

The default window size is 500bp. Efficiency can be improved by using a larger window size and step size. Shorter core-CREs can then be predicted within larger predicted regions. Core-CRE prediction in MOCCA uses motif occurrence predictions in order to identify high-scoring motif occurrence clusters with precise positions. With the following command, SVM-MOCCA is applied genome-wide using a window size of 3kb and step size of 1kb, and core-CREs are predicted based on clusters of positively classified motif occurrences within the windows.

```
mocca -auto:FASTA KahnPcG.fa -genome:FASTA dmel-all-chromosome-r5.57.fasta -motif:XML motifs2019expho.xml -motif:PWM UmassPGFE_PSSM_pho_SOLEXA_5_20200302.pssm 0 -predict:GFF predictions.gff -predict:Wig predictions.wig -predict:core -wSize 3000 -wStep 1000
```


### SVM-MOCCA &ndash; Alternative kernels
--------------------------

The models we have trained so far are linear. Alternatively, SVMs support non-linear classification using kernel functions. Non-linear kernels can enable the modelling of combinatorial feature enrichment. A quadratic model can be trained by adding the `-k:quadratic` argument (see `mocca --help` for more information on SVM kernel functions).

```
mocca -auto:FASTA KahnPcG.fa -genome:FASTA dmel-all-chromosome-r5.57.fasta -motif:XML motifs2019expho.xml -motif:PWM UmassPGFE_PSSM_pho_SOLEXA_5_20200302.pssm 0 -predict:GFF predictions.gff -predict:Wig predictions.wig -predict:core -wSize 3000 -wStep 1000 -k:quadratic
```


### SVM-MOCCA &ndash; Custom data
--------------------------

In preceding sections, we have used assisted learning. The user can also manually specify all training/validation/calibration data and classes to use. We will here train SVM-MOCCA using data of four classes: PREs, dummy PREs, dummy genomic and coding sequences. For PREs and coding sequences, we include files based on the Kahn *et al.* (2014) Polycomb target sites (`tPREsKahn2014.fa` and `vPREsKahn2014.fa` for independent training and validation sequences) and coding sequences from the FlyBase *D. melanogaster* genome annotation R5.57 (`tCDS.fa`).

Data classes are specified using the argument `-class VALUE LABEL`, where `VALUE` is an integer that identifies the class and is used with the machine learning algorithm, and `LABEL` designates the class as positive (`+`) or negative (`-`). Training data can be specified with FASTA files using the argument `-train:FASTA CLASS MODE`, where `CLASS` is a class value, and `MODE` is a training mode (`full` for the full sequences, or `windows` for training in windows). Training data can also be generated with `-train:MC PATH NUM LEN CLASS MODE ORDER`. See `mocca --help` for more information. Validation and calibration data is similarly specified with arguments with the `-train`-prefix replaced with `-validate` and `-calibrate` prefixes, and modes removed. We here use 11802 negative validation and calibrations sequences, which is the size of the genome divided by 3000 and scaled by the ratio of positives used for calibration to the full set (50/170).

When training without assisted learning, we additionally require setting the features to use for SVM-MOCCA. We will use local dinucleotide frequencies (`-f:MOCCA:DNT`) and local motif occurrence frequencies (`-f:MOCCA:nOcc`).

See `mocca --help` for more details regarding all other arguments.

```
mocca -motif:XML motifs2019.xml -wSize 3000 -wStep 1000 -C:SVM-MOCCA -f:MOCCA:DNT -f:MOCCA:nOcc -k:quadratic -class "PREs" 1 + -class "Dummy genomic" -1 - -class "Dummy PREs" -2 - -class "Coding sequences" -3 - -train:FASTA tPREsKahn2014.fa 1 full -train:MC dmel-all-chromosome-r5.57.fasta 110 3000 -1 full 4 -train:MC tPREsKahn2014.fa 110 3000 -2 full 4 -train:FASTA tCDS.fa -3 full -validate:FASTA vPREsKahn2014.fa + -validate:MC dmel-all-chromosome-r5.57.fasta 11802 3000 - 4 -calibrate:FASTA vPREsKahn2014.fa + -calibrate:MC dmel-all-chromosome-r5.57.fasta 11802 3000 - 4 -calibrate:precision 0.8 -genome:FASTA dmel-all-chromosome-r5.57.fasta -predict:GFF predictions.gff -predict:Wig predictions.wig -predict:core
```


### RF-MOCCA &ndash; Assisted training
--------------------------

In addition to including a polished implementation of SVM-MOCCA, MOCCA includes a derivative method based on Random Forests &ndash; RF-MOCCA. Assisted training of RF-MOCCA is achieved almost identically to with SVM-MOCCA, except that we add the argument `-C:RF-MOCCA`. We can also set the number of trees to use per motif classifier with the argument `-RF:trees COUNT`, where `COUNT` is the number of trees.

```
mocca -auto:FASTA KahnPcG.fa -genome:FASTA dmel-all-chromosome-r5.57.fasta -motif:XML motifs2019expho.xml -motif:PWM UmassPGFE_PSSM_pho_SOLEXA_5_20200302.pssm 0 -predict:GFF predictions.gff -predict:Wig predictions.wig -predict:core -wSize 3000 -wStep 1000 -C:RF-MOCCA -RF:trees 300
```


#### CPREdictor &ndash; Assisted training
--------------------------

Similarly, changing to training the CPREdictor model is as simple as adding the following argument: `-C:CPREdictor`. The PREdictor (Ringrose *et al.* 2003) is a binary classifier. The `-auto:FASTA` pipeline accordingly changes the negative training set generation to use one generated class of negatives &ndash; dummy-CREs.

```
mocca -auto:FASTA KahnPcG.fa -genome:FASTA dmel-all-chromosome-r5.57.fasta -motif:XML motifs2019.xml -predict:GFF predictions.gff -predict:Wig predictions.wig -C:CPREdictor
```


#### Models with user-specified feature spaces &ndash; Assisted training
--------------------------

MOCCA supports training log-odds, Support Vector Machine (SVM) and Random Forest (RF) models with a variety of motif-based feature sets. For a complete listing of feature sets implemented in MOCCA, see `mocca --help`.

To train a classic SVM, we use the argument `-C:SVM`. To use the motif occurrence pair frequency feature set with a pairing distance cutoff of 220 bp, we add the argument `-f:nPair 220`.
```
mocca -auto:FASTA KahnPcG.fa -genome:FASTA dmel-all-chromosome-r5.57.fasta -motif:XML motifs2019expho.xml -predict:GFF predictions.gff -predict:Wig predictions.wig -C:SVM -f:nPair 220
```

Simialrly, to instead train a classic RF, we use the argument `-C:RF`.
```
mocca -auto:FASTA KahnPcG.fa -genome:FASTA dmel-all-chromosome-r5.57.fasta -motif:XML motifs2019expho.xml -predict:GFF predictions.gff -predict:Wig predictions.wig -C:RF -f:nPair 220
```


#### Abbreviations
--------------------------
 * CRE: *Cis*-regulatory element
 * PRE: Polycomb/Trithorax Response Element
 * SVM: Support Vector Machine
 * RF: Random Forest
 * MOCCA: Motif Occurrence Combinatorics Classification Algorithm


-------------------------------------------------

#### References

 * Bredesen *et al.* 2019: https://academic.oup.com/nar/article/47/15/7781/5538007
 * Kahn *et al.* 2014: https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1004495
 * Gnocis - Bredesen *et al.* 2020 (?): Manuscript in preparation
 * FlyBase - : https://academic.oup.com/nar/article/41/D1/D751/1051942
 * Ringrose *et al.* 2003: https://www.sciencedirect.com/science/article/pii/S153458070300337X
 * Fly Factor Survey - https://academic.oup.com/nar/article/39/suppl_1/D111/2508103

