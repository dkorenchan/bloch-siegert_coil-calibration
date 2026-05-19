# TNMR_sequenceFiles

Tecmag TNMR data files (`.tnt`) and their associated pulse sequence files (`.tps`) for performing Bloch-Siegert coil calibration experiments on a low-field MR system.

In all sequences an off-resonance Bloch-Siegert pulse is applied during the echo period. The resulting signal phase encodes the local B₁ field strength via the Bloch-Siegert shift.

| File pair | Description |
|---|---|
| `BS1HprofileCal_ActiveTR` | Single spin-echo sequence for mapping the <sup>1</sup>H transmit coil profile as a function of offset frequency |
| `BS1HprofileCalMultiEcho_ActiveTR` | Rapid multi-echo variant of the <sup>1</sup>H coil profile calibration sequence |
| `BScalXcoil_ActiveTR` | Single spin-echo sequence for calibrating the X-nuclear transmit coil amplitude |
| `BScalXcoilMultiEcho_ActiveTR` | Rapid multi-echo variant of the X-nuclear coil calibration sequence |
