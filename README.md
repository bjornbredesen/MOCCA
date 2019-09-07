# MOCCA
Copyright Bjørn Bredesen, 2011-2019

-------------------------------------------------

## About
MOCCA (Motif Occurrence Combinatorics Classification Algorithms) is a suite for modelling DNA cis-regulatory element sequences.
Several types of motif-based models are included within MOCCA: a reimplementation of the PREdictor (Ringrose *et al.* 2003), implementations of the Dummy PREdictor and SVM-MOCCA (Bredesen *et al.* 2019), as various motif-based kernel functions that can be combined with log-odds and Support Vector Machine models.

#### References
 * Bredesen *et al.* 2019: https://academic.oup.com/nar/article/47/15/7781/5538007
 * Ringrose *et al.* 2003: https://www.sciencedirect.com/science/article/pii/S153458070300337X

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

## Included dependency licenses

### LibSVM

Copyright (c) 2000-2013 Chih-Chung Chang and Chih-Jen Lin
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions
are met:

1. Redistributions of source code must retain the above copyright
notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright
notice, this list of conditions and the following disclaimer in the
documentation and/or other materials provided with the distribution.

3. Neither name of copyright holders nor the names of its contributors
may be used to endorse or promote products derived from this software
without specific prior written permission.


THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE REGENTS OR
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

### RapidXML

Copyright (c) 2006, 2007 Marcin Kalicinski

Permission is hereby granted, free of charge, to any person obtaining a copy 
of this software and associated documentation files (the "Software"), to deal 
in the Software without restriction, including without limitation the rights 
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies 
of the Software, and to permit persons to whom the Software is furnished to do so, 
subject to the following conditions:

The above copyright notice and this permission notice shall be included in all 
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL 
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER 
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, 
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS 
IN THE SOFTWARE.

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

