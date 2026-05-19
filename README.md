# Bloch-Siegert Calibration of Low-field MR Coils
This repository contains Tecmag TNMR pulse sequence files for performing low-field MR coil calibration using Bloch-Siegert shifts, as well as processing code for the data.

## Citations
Please cite the following papers if you use the pulse sequences and/or code provided:

## Repository Contents

### `TNMR_sequenceFiles/`
Tecmag TNMR pulse sequence files (`.tnt`) and their associated parameter files (`.tps`) for performing Bloch-Siegert coil calibration experiments on a low-field MR system.

| File pair | Description |
|---|---|
| `BS1HprofileCal_ActiveTR` | Single spin-echo sequence for mapping the <sup>1</sup>H transmit coil profile as a function of offset frequency |
| `BS1HprofileCalMultiEcho_ActiveTR` | Rapid multi-echo variant of the <sup>1</sup>H coil profile calibration sequence |
| `BScalXcoil_ActiveTR` | Single spin-echo sequence for calibrating the X-nuclear transmit coil amplitude |
| `BScalXcoilMultiEcho_ActiveTR` | Rapid multi-echo variant of the X-nuclear coil calibration sequence |

In all sequences an off-resonance Bloch-Siegert pulse is applied during the echo period. The resulting signal phase encodes the local B₁ field strength via the Bloch-Siegert shift.

---

### `MATLAB_processingScripts/`
MATLAB functions for loading and processing Bloch-Siegert calibration data acquired with the sequences above.

| Script | Description |
|---|---|
| `BScoilProfileCal.m` | Processes <sup>1</sup>H coil profile data from the single spin-echo sequence. Returns the measured B₁ amplitude (Hz) and, optionally, corrected pulse amplitudes targeting a desired B₁, as a function of offset frequency. |
| `BScoilProfileCal_CPMG.m` | Same as above, adapted for the CPMG multi-echo acquisition. |
| `BSXcoilCal.m` | Processes X-nuclear coil calibration data from the single spin-echo sequence. Returns the measured B₁ amplitude as a function of pulse amplitude (linear scale factor). |
| `BSXcoilCal_CPMG.m` | Same as above, adapted for the CPMG multi-echo acquisition. |

#### `MATLAB_processingScripts/TNMRfileloadScripts/`
Helper functions for parsing Tecmag `.tnt` data files. Used internally by the calibration scripts above; these do not need to be called directly. See the [README](MATLAB_processingScripts/TNMRfileloadScripts/README.md) in that subfolder for details.
