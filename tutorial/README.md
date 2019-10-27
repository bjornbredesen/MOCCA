
![alt text](../markdown/mocca.png "")

# MOCCA
Copyright Bjørn Bredesen, 2011-2019

-------------------------------------------------

## Tutorial

### Assisted training
In this section, we will train SVM-MOCCA to predict Polycomb/Trithorax Response Elements (PREs) using Polycomb targets from the Kahn *et al.* (2014) experimental set. It is assumed that MOCCA has been installed. We extracted the Polycomb target regions determined by Kahn *et al.* (2014) from their Supplementary Table 3, we resized all regions to a length of 3kb each, and we extracted the underlying sequences from the *Drosophila melanogaster* genome release 5, using Cistem. The resulting sequences are supplied in the file `KahnPcG.fa`.

Sequence models in MOCCA make use of sequence motifs in order to learn distinguishing features of the different sequence classes. Accordingly, the motifs to be used must be specified. The simplest way of supplying MOCCA with a set of known motifs is by specifying them in IUPAC nucleotide codes (https://www.bioinformatics.org/sms/iupac.html) in an XML-file. We have supplied a sample XML-file, with the motifs that we used previously for genome-wide prediction of PREs (Bredesen *et al* 2019). The motifs are given in the file `motifs2019.xml`.

MOCCA uses binary or multi-class machine learning, and thus requires training data of at least two different classes. However, MOCCA supplies a mode that automatically constructs negative training data based on an input genome, given in FASTA format. We will train our first model using this mode. We will download the *Drosophila melanogaster* genome from FlyBase.

 * **Either**: Download and unpack the genome sequence to the folder with the tutorial data: ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r5.57_FB2014_03/fasta/dmel-all-chromosome-r5.57.fasta.gz
 * **Or**: In a Linux terminal, navigate to the folder with the tutorial data, and execute the following command: `wget ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r5.57_FB2014_03/fasta/dmel-all-chromosome-r5.57.fasta.gz && gunzip dmel-all-chromosome-r5.57.fasta.gz`

With this set up, we are ready to train SVM-MOCCA. In a terminal, navigate to the tutorial folder, and run:
```
mocca -auto:FASTA KahnPcG.fa -genome:FASTA dmel-all-chromosome-r5.57.fasta -motif:XML motifs2019.xml 
```

MOCCA will then do the following:
 * Load the input files.
 * Divide the training sequences (`KahnPcG.fa`) into independent training and test sets.
 * Calculate i.i.d. distributions for the genome and training sequences, and generate negative training and test data.
    - The negative test data is 100 times larger than the positive test data, in order to reflect the expected imbalance in the genome.
 * Generate threshold calibration data.
    - This is the size of the genome.
 * Train an SVM-MOCCA model, with three classes (CRM, dummy CRM and dummy genomic).
 * Calibrate the threshold for an expected genome-wide precision of 80% (can be changed).
 * Print out validation statistics.
 * *Note*: This may take some time to execute.

We usually also want to save the predictions. In order to do so, at least one additional argument must be specified. In order to output genome-wide predictions as non-overlapping regions in the General Feature Format (https://www.ensembl.org/info/website/upload/gff.html), add `-predict:GFF predictions.gff`. To output all genome-wide window scores, add `-predict:Wig predictions.wig`.
```
mocca -auto:FASTA KahnPcG.fa -genome:FASTA dmel-all-chromosome-r5.57.fasta -motif:XML motifs2019.xml -predict:GFF predictions.gff -predict:Wig predictions.wig
```

After running this, MOCCA should have generated the files `predictions.gff` and `predictions.wig` in the same folder.

The SVM-MOCCA model we here trained is a linear model. A quadratic model can be trained by adding the `-k:quadratic` argument (see `mocca --help` for more information on SVM kernel functions).


#### Remarks
--------------------------

 * The negative sequences generated with the `-auto:FASTA` option are generated with an i.i.d. sequence model. I.i.d. nucleotide models do not preserve motif occurrence information, and the generated sequences are thus likely to be too null to yield high precision. In Bredesen *et al.* (2019), we used 4th order Markov chains to generate negatives.

 * If prediction is to be performed on only a subset of the chromosomes specified in the input genome FASTA file, the FASTA file must be filtered in advance of prediction.


#### Assisted training of CPREdictor
--------------------------

Changing to training the CPREdictor model is as simple as adding the following argument: `-C:CPREdictor`. The PREdictor (Ringrose *et al.* 2003) is a binary classifier. The `-auto:FASTA` pipeline accordingly changes the negative training set generation to use one generated class of negatives.

```
mocca -auto:FASTA KahnPcG.fa -genome:FASTA dmel-all-chromosome-r5.57.fasta -motif:XML motifs2019.xml -predict:GFF predictions.gff -predict:Wig predictions.wig -C:CPREdictor
```


#### Assisted training of models with user-specified feature spaces
--------------------------

MOCCA supports training log-odds or Support Vector Machine (SVM) models with a variety of motif-based feature sets. For a complete listing of feature sets implemented in MOCCA, see `mocca --help`.

A Support Vector Machine with motif pair occurrence frequencies as features can be trained using the following command:
```
mocca -auto:FASTA KahnPcG.fa -genome:FASTA dmel-all-chromosome-r5.57.fasta -motif:XML motifs2019.xml -predict:GFF predictions.gff -predict:Wig predictions.wig -C:SVM -f:nPair 220
```


-------------------------------------------------

#### References

 * Bredesen *et al.* 2019: https://academic.oup.com/nar/article/47/15/7781/5538007
 * Kahn *et al.* 2014: https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1004495
 * Cistem - Bredesen *et al.* 2020 (?): Article to be submitted to journal
 * FlyBase - : https://academic.oup.com/nar/article/41/D1/D751/1051942
 * Ringrose *et al.* 2003: https://www.sciencedirect.com/science/article/pii/S153458070300337X

