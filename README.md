# Extended Multi-Particle Correlation (EMPC) via HIJING model
A C++ and ROOT-based framework integrated with the HIJING: Heavy Ion Jet INteraction Generator model to efficiently calculate $N$-particle azimuthal correlations and flow cumulants in heavy-ion collisions.

## Overview
This repository contains a modified C++ implementation of the Extended Multi-Particle Correlation (EMPC) technique, integrated directly into the HIJING model. The primary codebase has been translated and structured in C++ using the ROOT framework for robust data handling, replacing traditional standalone Fortran analysis workflows.

## General Functionality
### 1. Generating Q-Vectors and Multiparticle Correlators
The core mathematical engine is driven by the `GeneralQCorrelator` and `WrapperQCorrEvent` classes. Instead of iterating through nested loops to calculate multiparticle correlations, the framework uses a generalized $Q$-vector approach. 
* Tracks are processed to calculate phase factors (e.g., $e^{i n \phi}$) which are temporarily cached for the event.
* The algorithm dynamically constructs block sums and evaluates mathematical partitions to calculate exact integrated and differential multi-particle cumulants up to 6-particle correlations ($d_{n_1, n_2, \dots}$).
* It seamlessly supports both single-event and two-subevent (with $\eta$ gaps) correlation methods to suppress non-flow effects.

### 2. Filling Histograms
Output and data storage are managed by the `Wrapper` and `HistoAutoFile` classes. 
* Hundreds of `TProfile` histograms are automatically booked in memory to track harmonic correlations across varying centralities (based on participant counts) and transverse momentum ($p_T$) bins.
* As the HIJING model generates events, kinematic properties (like $p_T$, $\eta$, $\phi$) are extracted and fed into the `WrapperQCorrEvent` correlator.
* After the $Q$-vectors and sub-event combinations are evaluated, the resulting real components and event weights are filled directly into the corresponding `TProfile` bins. 

### 3. Combining Fortran and C++ Codes
This framework utilizes `CMake` to bind the original Fortran HIJING generators with the modern C++ analysis backbone. 
* Fortran subroutines are declared as `extern "C"` blocks in the main C++ file, allowing direct memory access to Fortran common blocks.
* The CMake configuration (`CMakeLists.txt`) natively handles the cross-language compilation, linking the `gfortran` compiler outputs with C++17 and ROOT libraries to produce a single unified executable (`HIJrun`).

## Citation Request
If you utilize this modified C++ framework or the integrated EMPC correlation algorithms in your research, please ensure you cite the original author and the foundational paper:
> *Niseem Magdy, 
> *[e-Print: 2405.19169]*

> *Niseem Magdy,
> *[e-Print: 2607.17449]*

*Please also link back to this repository in your project documentation.*
