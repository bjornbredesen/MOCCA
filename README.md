# MOCCA
Copyright Bjørn Bredesen, 2011-2019

## About
MOCCA is a suite for modelling DNA cis-regulatory element sequences.

## Features
 - SVM-MOCCA (the Support Vector Machine Motif Occurrence Combinatorics Classification Algorithm)
 - CPREdictor
 - Dummy PREdictor
 - IUPAC motif occurrence parsing Finite State Machine

## Using
MOCCA accepts the following arguments:
 - `-C:SVM-MOCCA`
   * Use SVM-MOCCA.
 - `-C:SVM-MOCCA:C-SVC`
   * Use SVM-MOCCA with $C$-SVC.
 - `-C:SVM-MOCCA:nu-SVC`
   * Use SVM-MOCCA with $\nu$-SVC.
 - `-C:CPREdictor`
   * Use CPREdictor.
 - `-C:DummyPREdictor`
   * Use dummy PREdictor.
 - `-k:linear`
   * Use linear kernel.
 - `-k:quadratic`
   * Use quadratic kernel.
 - `-k:cubic`
   * Use cubic kernel.
 - `-k:RBF`
   * Use radial basis function kernel.
 - `-SVM:gamma VALUE`
   * Set the SVM $\gamma$ kernel parameter.
 - `-SVM:c0 VALUE`
   * Set the SVM $c_0$ kernel parameter.
 - `-SVM:p VALUE`
   * Set the SVM $\rho$-hyperparameter.
 - `-SVM:C VALUE`
   * Set the SVM $C$-hyperparameter.
 - `-SVM:nu VALUE`
   * Set the SVM $\nu$-hyperparameter.
 - `-wSize SIZE`
   * Sets the window size to `SIZE`.
 - `-wStep SIZE`
   * Sets the window step size to `SIZE`.
 - `-motif:kmer K`
   * Adds the set of all motifs of length `K` to the set of motifs.
 - `-motif:Random N LEN`
   * Adds `N` randomly generated motifs of length `LEN`.
 - `-motif:IUPAC NAME MOTIF MISMATCHES`
   * Adds an IUPAC-motif with the name `NAME`, and sequence `MOTIF`, with `MISMATCHES` mismatches allowed.
 - `-motif:XML PATH`
   * Adds motifs specified in an XML-file.
 - `-f:MOCCA:nOcc`
   * Adds motif occurrence frequency features to SVM-MOCCA.
 - `-f:MOCCA:DNT`
   * Adds dinucleotide features to SVM-MOCCA.
 - `-f:MOCCA:GC`
   * Adds GC-content features to SVM-MOCCA.
 - `-train:FASTA PATH CLASS MODE`
   * Adds a training sequence file.
   * `PATH`: Path to FASTA file.
   * `CLASS`: A class ID, defined with `-class`, or one of the pre-specified binary classes: `+` for positive or `-` for negative.
   * `MODE`: Can be `win`, for training with all windows within each training sequence file, or `full`, for training with the full sequences.
 - `-validate:FASTA PATH CLASS`
   * Adds a validation sequence file.
   * `PATH`: Path to FASTA file.
   * `CLASS`: A class ID, defined with `-class`, or one of the pre-specified binary classes: `+` for positive or `-` for negative.
 - `-validate:outSCTable PATH`
   * Writes scores and classes for each validation sequence to a tab-separated table file, `PATH`.
 - `-class NAME VALUE FLAG`
   * Registers a sequence class with name `NAME`, value/ID `VALUE`, and flag `FLAG` (\"+\" for positive, or \"-\" for negative)

### Example usage
`./bin/mocca -motif:XML ./data/motifs.2003.xml -motif:IUPAC Combgap GTGT 0 -f:MOCCA:nOcc -f:MOCCA:DNT -train:FASTA ./data/T2003_PRE.fasta + full -train:FASTA ./data/T2003_NonPRE.fasta - full`

