# Compendium of CIN repository hub

## Structure and content of the five repositories associated with this study

This repository works as the central code hub for "A pan-cancer compendium of chromosomal instability", Drews et al. (2022). It contains all code necessary to reproduce the figures (main, extended and supplement). Linked to this repository are four other repositories which contain methods that were used during this study, but they can also be adapted to be used by others on their own data. These are:
  <!---* [Signature discovery](https://github.com/markowetzlab/SignatureDiscovery)-->

  * [Signature discovery](https://github.com/markowetzlab/CINSignatureDiscovery)

   This repository contains code to identify pan-cancer signatures of CIN. The approach takes ASCAT copy number profiles as input and generates the final signature definitions and activities.

  <!---* * [Signature quantification](https://github.com/markowetzlab/SignatureQuantification)-->

  * [Signature quantification](https://github.com/markowetzlab/CINSignatureQuantification)

   This repository contains code to quantify CIN signatures given an absolute copy number profile.

  <!--- * [Signature interpretation](https://github.com/macintyrelab/SignatureInterpretation) -->

  * [Signature interpretation](https://github.com/macintyrelab/CINSignatureInterpretationMatrix)

   This repository contains code to generate the signature interpretation matrix which is used to assist in determining signature aetiologies.

  <!--- * [Genome simulation](https://github.com/markowetzlab/CINGenomeSimulation) -->

  * [Genome simulation](https://github.com/markowetzlab/CINSignatureGenomeSimulation)

   This repository contains code to simulate CIN producing mechanism by generating copy number profiles. From these profiles, signatures can be identified and compared to the expected signatures.



## Linking code to sections in the manuscript

The code in this repository largely mirrors the structure of the methods section of the manuscript. The following table details the sections of the Methods and where to find the associated code.

| #   | Section                               | Details     | Location            |
| --- | ------------------------------------- | ----------- | --------------------|
| 1   | Generation of copy number profiles    | Contains sample profiles | [Vanloolab/ASCAT](https://github.com/VanLoo-lab/ascat/tree/master/ReleasedData/TCGA_SNP6_hg19) |
| 2   | Copy number feature encoding          | .           | [markowetzlab/CINSignatureDiscovery](https://github.com/markowetzlab/CINSignatureDiscovery) |
| 3   | Deriving signatures                   | .           | [markowetzlab/CINSignatureDiscovery](https://github.com/markowetzlab/CINSignatureDiscovery) |
| 4   | Assessing feature encoding performance| Simulating CIN | [markowetzlab/CINGenomeSimulation](https://github.com/markowetzlab/CINSignatureGenomeSimulation) |              |
| 5   | Robustness analysis                   | Activities & Definitions | This repo |
| 5.1 | Signature activity stability assessment via simulation | .  | This repo |
| 5.2 | Signature stability across genomic technologies | .         | This repo and [markowetzlab/CINSignatureSWGSBinsize](https://github.com/markowetzlab/CINSignatureSWGSBinsize) |
| 5.3 | Signature definition stability assessment via introducing noise | . | This repo |
| 6   | Enabling signature interpretation     | Shrinkage  | [macintyrelab/CINSignatureInterpretationAnalyses](https://github.com/macintyrelab/CINSignatureInterpretationAnalysis)|
| 6.1 | Feature shrinkage is a general property of NMF | .| [macintyrelab/CINSignatureInterpretationAnalyses](https://github.com/macintyrelab/CINSignatureInterpretationAnalysis)|
| 6.2 | Validation of interpretation matrix using simulations| .| [macintyrelab/CINSignatureInterpretationAnalyses](https://github.com/macintyrelab/CINSignatureInterpretationAnalysis) |
| 6.3 | Procedure to create signature interpretation matrix | Interpreting CNSigs | [macintyrelab/CINSignatureInterpretationMatrix](https://github.com/macintyrelab/CINSignatureInterpretationMatrix)|
| 7   | Identification of putative signature aetiologies | . | This repo|
| 7.1 | Linking signature activities to mutated genes  | . | This repo  |
| 7.2 | Heatmap activity by cancer type                | . | This repo  |
| 7.3 | IHR in CX2, CX3 and CX5                        | . | This repo  |
| 7.4 | Estimating the number of CNAs produced by a signature| . | This repo |
| 8   | Signatures for drug response prediction and target identification | .| [macintyrelab/CINSignatureBiomarkerAnalysis](https://github.com/macintyrelab/CINSignatureBiomarkerAnalysis)|
| 8.1 | Calculating signature activities for cell lines     | .| [markowetzlab/CINSignatureCellLines](https://github.com/markowetzlab/CINSignatureCellLines)|
| 9   | Building a clin. classifier bases on sig. activities| .| This repo |
| 10  | Supplementary figures                               | .| This repo |
| 11  | Miscellaneous (analyses without figures)            | .| This repo |
| 12  | Testing single-count feature encoding               | .| This repo |



## Linking code to manuscript figures (main, extended and supplement)

The directories above outline where to find the code associated with each Methods section. The table below explicitly lists which code reproduces each figure in the manuscript:

| Figure type | Figure number  | Details        | Code location      |
| ----------- | -------------- | -------------- | ------------- |
| Main        | 1              | No code needed | NA            |
| Main        | 2              | .              | This repository, directory Section 7 -> Section 7.3 |             |
| Main        | 3              | .              | [macintyrelab/CINSignatureBiomarkerAnalysis](https://github.com/macintyrelab/CINSignatureBiomarkerAnalysis)|
| Main        | 4              | .              | This repository, directory Section 9     |
| Extended    | 1              | .              | This repository, directory Section 1-4 Signature Discovery |
| Extended    | 2              | No code needed | NA            |
| Extended    | 3              | No code needed | NA            |
| Extended    | 4              | .              | [macintyrelab/CINSignatureInterpretationMatrix](https://github.com/macintyrelab/CINSignatureInterpretationMatrix)|
| Extended    | 5              | .              | This repository, directory Section 5 -> Section 5.1 |
| Extended    | 6              | .              | This repository, directory Section 5 -> Section 5.2 |
| Extended    | 7              | No code needed | NA            |
| Extended    | 8              | .              | This repository, directory Section 7 -> Section 7.1 |
| Extended    | 9              | .              | This repository, directory Section 7 -> Section 7.4 |
| Extended    | 10             | .              | This repository, directory Section 9     |
| Supplement  | 1-32           | .              | This repository, directory Section 10    |
| Supplement  | 33-34          | .              | This repository, directory Section 5 -> Section 5.3 |
| Supplement  | 40             | No code needed | NA  |
| Supplement  | 41             | .              | This repository, directory Section 10    |
| Supplement  | 42             | .              | This repository, directory Section 5.3 -> Noisy signatures |
| Supplement  | 43 & 51        | .              | [macintyrelab/CINSignatureInterpretationAnalysis](https://github.com/macintyrelab/CINSignatureInterpretationAnalysis)|
| Supplement  | 44 - 47 & 52   | Genome Simulation | This repository, directory Section 10   |
| Supplement  | 48             | . | This repository, directory Section 10   |
| Supplement  | 49 & 50        | .              | This repository, directory Section 12 |
| Supplement  | 54             | No code needed | NA  |

NB:

  * Section 11 contains scripts covering preparatory and additional analyses not captured by any figure.

## Contact

If you experience any issues or have questions about the code, please open a Github issue with a minimum reproducible example. For questions around collaborations or sensitive patient data, please contact us directly at Florian Markowetz <Florian.Markowetz@cruk.cam.ac.uk> and Geoff Macintyre <gmacintyre@cnio.es>.

## Licence
The contents of this repository are copyright (c) 2022, University of Cambridge and Spanish National Cancer Research Centre (CNIO).

The contents of this repository are published and distributed under the GAP Available Source License v1.0 (ASL). 

The contents of this repository are distributed in the hope that it will be useful for non-commercial academic research, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the ASL for more details. 

The methods implemented in the code are the subject of pending patent application GB 2114203.9.

Any commercial use of this code is prohibited.
