# ZtoEe-Analysis

In this study, a C language implementation of CompHep infrastructure was used to analyse the ZtoEe process. Kinematic distributions of electron positron and Z boson were studied. (pp->Z->eE). First of all, the files (txt file) produced in CompHep were converted to "lhe" format. This "lhe" file, which was also created with the ExRootLHEFConverter command, was converted to "ROOT" format.

Using the command LHEF->MakeClass("ZtoEeAnalysis"), the C and Header file of the ROOT file were created. To perform the analysis, the following commands are applied after the desired analysis codes are written in the C file.

.L ZtoEeAnalysis.C

ZtoEeAnalysis t

t.Loop()

Afterwards, the analysis ROOT files for electrons, positrons, and Z bosons are ready. The Breit-Wigner Function is also applied to the Z boson.

## Mass Paiz Z Histogram

![Screenshot_20220110_165821](https://user-images.githubusercontent.com/62266472/206757851-e3cad68a-2099-4dc8-8606-b6ad09c6a133.png)

