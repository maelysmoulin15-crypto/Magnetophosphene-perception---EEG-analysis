EEG Signals — Magnetic Field (MF) Exposure
1. Slicing Raw EEG Signals into Epochs
This repository is inspired by an internal MATLAB script by Julien Modolo (circa 2014) for MF–EEG experiments, which detects exposure onsets on a clean reference channel and applies the same epochs to all EEG channels in the exact randomized order.
The present code is an independent re-implementation focused on raw, robust epoching (EEGLAB .cnt I/O, auto channel selection, explicit sham placement, NaN padding policy, subject-wise overrides), with no filtering nor spectral analysis performed here.

2. EEG DSP (30–80 Hz) with Notches & 1/f Correction
What this script does
- Loads per-epoch EEG files.
- Applies your custom 30–80 Hz band-pass (filtrealex).
- Applies stable FIR notches at 40/50/60 Hz (linear-phase, deep attenuation).
- Computes the power spectrum (FieldTrip ft_freqanalysis) and FOOOF aperiodic (1/f) component (output='fooof_aperiodic').
- Derives oscillatory power = total − aperiodic, then computes gamma-band metrics (absolute power and power-weighted mean frequency) while excluding notched bins.
- Exports a results table (CSV + MAT) and summary PSD figures.

Requirements
- FieldTrip (recent version with FOOOF wrapper; addpath + ft_defaults).
- The custom function filtrealex available on your MATLAB path (3 files).
- Epoch files in .mat format containing at least:
- EEG_segment (vector), and either --> File naming pattern: Sujet%02d_%s_%dHz_%s_%dmT_trial%d.mat (subject, coil tag, stim frequency, electrode label, intensity mT, trial index)
- Optional: a behavioral/press table (CSV/XLSX) with columns named like Sujet/Subject, intensite/Intensity, essai/Trial, press/Press/Pression.

Configure paths : Edit the first section:
- addpath('/Users/maelys/Documents/MATLAB/') → folder that contains filtrealex
- addpath('/Users/maelys/Documents/MATLAB/fieldtrip') → FieldTrip root
- epoch_dir, results_dir, press_xlsx as needed
