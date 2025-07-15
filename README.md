# Telomere Analysis in Intestinal Crypts (Specific ROIs) of Liver Tissue

This repository contains a 3D image analysis pipeline for the quantification of telomeres in intestinal crypts of liver tissue using confocal microscopy data. The pipeline includes preprocessing, segmentation, and quantitative analysis at both the nuclear and telomere object levels.

## Overview

The pipeline performs the following steps:

1. **Preprocessing**: Applies manually curated ROIs to z-stacks and removes background outside crypt regions.
2. **Segmentation**: Uses custom-trained 3D deep learning models in Cellpose to segment nuclei (DAPI) and telomeres (TRF1).
3. **Quantification**:
   - **Per nucleus**: Computes number of telomeres, mean/sum intensity, and stratifies foci by intensity percentiles.
   - **Per telomere**: Measures intensity, volume, and spatial association with nuclei.

## Requirements

- Fiji (ImageJ) with Groovy scripting support
- MCIB3D library
- Cellpose (for training and inference)
- Java 8+

## Usage

1. Run `preprocess.groovy` to apply ROIs and clean background.
2. Segment nuclei and telomeres using Cellpose with your trained models.
3. Run `3DIntensitycPerNucleus.groovy` and `3DIntensitycPerObject.groovy` for quantification.
4. Results are saved in `output/csv/`.


