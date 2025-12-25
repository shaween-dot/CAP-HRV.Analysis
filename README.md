# CAP-HRV.Analysis
MATLAB-based framework for integrated EEG–ECG analysis of CAP and HRV during sleep.

This repository provides a MATLAB-based, modular EEG–ECG processing framework for:
-CAP detection during slow-wave sleep (SWS)
-Epoch-level ECG segmentation
-HRV feature extraction
-Integrated analysis of sleep microstructure and autonomic dynamics
-The code is intended as a technology and methods framework, not as a novel HRV algorithm implementation.

HRV features in this pipeline are computed using routines derived from the open-source PhysioZoo platform:
Behar JA et al., PhysioZoo: a novel open-access platform for heart rate variability analysis, Frontiers in Physiology, 2018.
The PhysioZoo algorithms themselves are not re-invented here. Instead, they are embedded, adapted, and integrated into an EEG–ECG workflow that enables CAP-aligned HRV analysis.

Validation included comparison of selected HRV metrics between:
-HRV computed inside the integrated MATLAB pipeline
-HRV computed using the PhysioZoo stand-alone platform

This framework provides:
-A unified EEG–ECG processing pipeline
-Automated CAP-aligned HRV analysis
-Reproducible, modular, open MATLAB code suitable for research use

Intended Use

This code is intended for:
-Sleep physiology research
-Methodological development
-Multimodal EEG–ECG analysis
-Reproducible feature extraction pipelines
