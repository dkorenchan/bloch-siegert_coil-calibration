# MATLAB_processingScripts

MATLAB functions for loading and processing Bloch-Siegert calibration data acquired with the sequences in `TNMR_sequenceFiles/`.

| Script | Input sequence | Description |
|---|---|---|
| `BScoilProfileCal.m` | `BS1HprofileCal_ActiveTR` | Processes <sup>1</sup>H coil profile data from the single spin-echo sequence. Returns the measured B₁ amplitude (Hz) and, optionally, corrected pulse amplitudes targeting a desired B₁, as a function of offset frequency. |
| `BScoilProfileCal_CPMG.m` | `BS1HprofileCalMultiEcho_ActiveTR` | Same as above, adapted for the CPMG multi-echo acquisition. |
| `BSXcoilCal.m` | `BScalXcoil_ActiveTR` | Processes X-nuclear coil calibration data from the single spin-echo sequence. Returns the measured B₁ amplitude as a function of pulse amplitude (linear scale factor). |
| `BSXcoilCal_CPMG.m` | `BScalXcoilMultiEcho_ActiveTR` | Same as above, adapted for the CPMG multi-echo acquisition. |

## `TNMRfileloadScripts/`
Helper functions for parsing Tecmag `.tnt` data files. Used internally by the calibration scripts above; these do not need to be called directly. See the [README](TNMRfileloadScripts/README.md) in that subfolder for details.
