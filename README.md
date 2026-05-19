# Bloch-Siegert Calibration of Low-field MR Coils
This repository contains Tecmag TNMR pulse sequence and data files for performing low-field MR coil calibration using Bloch-Siegert shifts, as well as processing code for the data.

> **Note:** The contents of this repository were originally designed for use on a 6.5 mT low-field MR system.

## Citations
Please cite the following papers if you use the pulse sequences and/or code provided:

## Repository Contents

### `TNMR_sequenceFiles/`
Tecmag TNMR data files (`.tnt`) and their associated pulse sequence files (`.tps`) for performing Bloch-Siegert coil calibration experiments on a low-field MR system. The provided pulse sequence files can be used for two main functions:

1. Calibration of an X-nuclear coil in a dual-tuned <sup>1</sup>H/X probe, using Bloch-Siegert shifts detected on <sup>1</sup>H
2. Off-resonance transmit B₁ profile mapping of a <sup>1</sup>H coil, using Bloch-Siegert shifts

See the [README](TNMR_sequenceFiles/README.md) in that subfolder for details.

### `MATLAB_processingScripts/`
MATLAB functions for loading and processing Bloch-Siegert calibration data acquired with the sequences above, including helper functions for parsing Tecmag `.tnt` data files. See the [README](MATLAB_processingScripts/README.md) in that subfolder for details.
