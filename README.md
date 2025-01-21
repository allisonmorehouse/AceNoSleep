# AceNoSleep

## Overview

This script processes and analyzes EEG data by integrating information from multiple file types, including sleep EDF files and Kubios HRV outputs. The analysis includes filtering EEG signals, detecting HR bursts, and calculating power in various EEG frequency bands.

---

## Requirements

### Input Files

1. **EDF file**: Contains EEG data.
2. **MAT file**: Kubios HRV output containing R-peak data from ECG analysis.

### Software

- MATLAB
- Signal Processing Toolbox

---

## Instructions

1. **File Preparation**
    
    Ensure the following:
    
    - Files are named consistently so they can be matched automatically (e.g., by a subject identifier).
    - Input files are placed in their respective directories:
        - `edfFolder`: Directory for EDF files.
        - `hrvFolder`: Directory for Kubios HRV outputs.
2. **Parameter Configuration**
    
    Adjust the parameters as needed:
    
    - `segmin`: Minimum duration (in minutes) of a stable stage for inclusion.
    - `hrbwin`: Window size (in seconds) around detected HR bursts.
    - `avbin`: Bin size for averaging EEG powers.
    - Frequency bands for EEG power calculations (`swa`, `alpha`, `theta`, `sigma`) and HRV (`lf`, `hf`).
        - Delta (0.5 - 4 Hz) or swa
        - SlowDelta (0.5 - 1 Hz) or swa1
        - FastDelta (1 to 4 Hz) or swa2
        - Alpha (8-13 Hz)
        - Sigma (12-16 Hz)
        - Theta (4-7 Hz)
3. **Update Code Sections**
    
    Use the search function to locate and modify occurrences of `"update here"`:
    
    - **File paths**: Specify directories for input files and output results.
    - **Channel selection**: Update the EEG channels if desired channels are different from the following (`{'F3', 'F4', 'C3', 'C4', 'O1', 'O2', 'Cz'}`).
    - **File naming conventions**: Modify parsing logic if filenames deviate from the expected format.
4. **Run the Script**
    
    Execute the script in MATLAB. The script will:
    
    - Match and load corresponding files for each subject.
    - Process EEG and HRV data, including filtering and feature extraction.
5. **Output Files**
    
    Results will be saved in the directory specified by `outputFolder`. There will be one file per subject.
    

---

## Output Variables

The data is saved in a mat structure called “ace.”

Bins:

- bin 1 = -10 to -5 seconds before the HRB
- bin 2 = -5 to 0 seconds before the HRB
- bin 3 = 0 to 5 seconds after the HRB
- bin 4 = 5 to 10 seconds after the HRB

Variables:

- bin{Frequency}_{electrode}: Binned power data for different frequency bands and sleep stages
    - Each variable will have 4 columns (or bins) and each row is a heart rate burst
- bin{HRV variable}: Binned power data for HRV variable.
    - Each variable will have 4 columns (or bins) and each row is a heart rate burst
- av{Frequency}_{electrode}: Binned power data for different frequency bands.
    - Each variable will have 4 columns (or bins) and 1 row, which is the average of all of the heart rate burst
- HRB_{Frequency}_allCh: electrode x 1 cell
    - For each electrode - Power in the channel (row) x timepoints (column)
- bl{Frequency}_{electrode}: Baseline power for different frequency bands for the electrode specified
    - Result is one value
- HRB_idx = Indices of heart rate bursts for each sleep stage
