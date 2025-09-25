EEG Signals — Magnetic Field (MF) Exposure
Slicing Raw EEG Signals into Epochs

This repository is inspired by an internal MATLAB script by Julien Modolo (circa 2014) for MF–EEG experiments, which detects exposure onsets on a clean reference channel and applies the same epochs to all EEG channels in the exact randomized order.
The present code is an independent re-implementation focused on raw, robust epoching (EEGLAB .cnt I/O, auto channel selection, explicit sham placement, NaN padding policy, subject-wise overrides), with no filtering nor spectral analysis performed here.
