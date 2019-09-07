# MOCCA
Copyright Bjørn Bredesen, 2011-2019

-------------------------------------------------

## About
MOCCA (Motif Occurrence Combinatorics Classification Algorithms) is a suite for modelling DNA cis-regulatory element sequences.
Several types of motif-based models are included within MOCCA: a reimplementation of the PREdictor (Ringrose et al. 2003), implementations of the Dummy PREdictor and SVM-MOCCA (Bredesen et al. 2019), as various motif-based kernel functions that can be combined with log-odds and Support Vector Machine models.

#### References
Bredesen *et al* 2019: https://academic.oup.com/nar/article/47/15/7781/5538007
Ringrose *et al* 2003: https://www.sciencedirect.com/science/article/pii/S153458070300337X

-------------------------------------------------

## License

MIT License

Copyright (c) 2019 Bjørn André Bredesen

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

-------------------------------------------------

## Features
 - Models
     * Dummy PREdictor
     * CPREdictor
     * SVM-MOCCA (the Support Vector Machine Motif Occurrence Combinatorics Classification Algorithm)
     * Log-odds models with motif-based feature spaces
     * Support Vector Machines with motif-based feature spaces
 - Motif handling
     * Command-line specification of IUPAC motifs
     * Loading of IUPAC motifs from XML
     * Generation of random IUPAC motifs.
     * Full *k*-mer sets.
     * IUPAC motif occurrence parsing Finite State Machine
 - Feature spaces
     * Motif occurrence frequency spectrum
     * Motif pair occurrence frequency spectrum, with distance cutoff, and multiple distancing and overlap modes
     * Motif distancing kernels
     * Periodic motif occurrence kernels
     * Motif pairing kernels that incorporate positional information
 - Core usage features
     * Training with FASTA sequence files
     * Validation with FASTA sequence files
     * Saving sequence scores to table
     * Scoring of sequence files to Wiggle curves

### Example usage
`./bin/mocca -motif:XML ./data/motifs.2003.xml -motif:IUPAC Combgap GTGT 0 -f:MOCCA:nOcc -f:MOCCA:DNT -train:FASTA ./data/T2003_PRE.fasta + full -train:FASTA ./data/T2003_NonPRE.fasta - full`

