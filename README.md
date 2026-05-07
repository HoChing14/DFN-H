DFN-H

DFN-H is a Fortran 90/95 toolbox developed for efficient transient coupled flow and heat transport forward modeling in fractured rocks based on the TOUGHREACT simulator. It employs an analytical solution to solve heat transfer between fracture and the low-permeability rock matrix while neglecting fluid (water) flow between them. This eliminates the need to explicitly discretize the rock matrix, while retaining the essential thermo-hydraulic coupling processes that governs fractured reservoir behavior.

The repository contains DFN-H core solver (including the developed DFN-H source code files, a Python script for compling, and a file contains compling sequence of the source code files), four reproducible synthetic examples and their initialization and visualization codes, and reference materials for the published methods.

1 Scope

DFN-H currently supports:
Coupled flow and heat transport forward modeling in 2D and 3D fractured media with negligible rock matrix permeability in time-independent continuous injection conditions
COMSOL Multiphysics (DFM) forward modeling results preprocessing and visulization
Spatial and temporal fracture temperature distribution visualization

2 Repository layouts

tough_codes/ : DFN-H core code files, gfortran compling script (for Linux), and a compling dependence file (sequence.txt)
ana_solu/ : DFN-H accuracy verification example (single fracture analytical solution)
nume_1/ : DFN-H genericity verification example 1 (single fracture)
nume_2/ : DFN-H genericity verification example 2 (two orthogonal fractures)
nume_3/ : DFN-H efficiency verification example (stochastic fracture network)

3 Requirements

Linux
gfortran for Fortran 90/95 code compliation
Python (version 3.7+, mandatory packages: numpy, scipy, math, pandas, re, copy, matplotlib, os, subprocess, shutil, time)
Paraview

4 Quick start
4.1 clone codes:
git clone https://github.com/HoChing14/DFN-H.git
4.2 DFN-H core code compliation:
cd ../tough_codes
Before compliation, open the file "sequence.txt" in this directory and modify the first line to your username, then:
python main.py
The compliation succeeds when no information except "Compilation successful!" appears on the screen.

4.3 Run accuracy verification example
cd ../ana_solu
python main.py
If the run succeeds, output files named "spatial.txt" and "temporal.txt" (spatial and temporal temperature variation results of analytical solution and DFN-H numerical solution) are created or overwritten. The above results can be imported into Origin or Excel for visulization.

4.4 Run genericity and efficiency verification examples
The following commands executes the forward modeling of the corresponding example:
cd ../nume_1 (or nume_2, nume_3)
python main.py
If the run succeeds, output file named "cal_T.txt" (calculated temporal temperature variation) is created or overwritten, and can be imported into Origin or Excel for visulization.
