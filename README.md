# ECG-AF-Detection-EEG-Sleep-Staging

This repository contains a project focused on the analysis and classification of ECG and EEG signals for detecting atrial fibrillation (AF) and classifying sleep stages (Wake, S2, S3).

The project is divided into two main parts:
1. **ECG-based Detection of Atrial Fibrillation (AF)**
2. **EEG-based Sleep Staging**

---

## Table of Contents
1. [Introduction](#introduction)
2. [ECG-based Atrial Fibrillation Detection](#ecg-based-atrial-fibrillation-detection)
   - [Preprocessing](#preprocessing)
   - [Artifact Removal](#artifact-removal)
   - [Feature Extraction](#feature-extraction)
   - [Unsupervised Classification](#unsupervised-classification)
   - [Evaluation Metrics](#evaluation-metrics)
3. [EEG-based Sleep Staging](#eeg-based-sleep-staging)
   - [Segmentation](#segmentation)
   - [K-complex Detection](#k-complex-detection)
   - [Power Spectral Density Analysis](#power-spectral-density-analysis)
   - [Classification Models](#classification-models)
4. [Results](#results)
5. [Conclusion](#conclusion)

---

## Introduction

This project explores the use of signal processing and machine learning techniques for:
1. Detecting atrial fibrillation (AF) from ECG signals.
2. Classifying sleep stages (Wake, S2, S3) using EEG data.

The project employs preprocessing, feature extraction, and classification techniques for each task.

---

## ECG-based Atrial Fibrillation Detection

### Preprocessing
The preprocessing steps clean the raw ECG signals by:
- Time-domain and frequency-domain analysis to identify noise and artifacts.
- Use of filters to remove powerline interference, motion artifacts, and baseline wander.

![Add image of artifact removal results here]()

### Feature Extraction
We extracted features from the ECG signal, including:
- **Morphological features**: Signal energy, wavelet subband energies.
- **Heart-rate variability (HRV) features**: RMSSD, pRR50, and minRR.

### Unsupervised Classification
Using Principal Component Analysis (PCA) and K-means clustering, the ECG segments were classified into 'normal' and 'AF' categories.

### Evaluation Metrics
The performance of the classification was evaluated based on:
- **Accuracy**: 70.89%
- **Sensitivity**: 71.57%
- **Specificity**: 62.50%

---

## EEG-based Sleep Staging

### Segmentation
The EEG signal was segmented into 30-second intervals for further analysis.

![Add image of K-complex detection here]()

### Power Spectral Density Analysis
The Power Spectral Density (PSD) of the EEG signal was computed for 10 different frequency bands to analyze the signal's frequency characteristics during different sleep stages.

### Classification Models
A neural network classifier with two hidden layers was implemented to classify the EEG signal into sleep stages (Wake, S2, S3). Feature selection was performed using both manual methods and Linear Discriminant Analysis (LDA).

---

## Results

- **ECG Classification Results**: The unsupervised model achieved an accuracy of 70.89% for AF detection.
- **EEG Classification Results**: The neural network classifier achieved an accuracy of 72.79% with the full feature set and 81.80% with a reduced set of features.

---

## Conclusion

This project demonstrates the effectiveness of signal processing techniques combined with machine learning models for biomedical data analysis, specifically for detecting atrial fibrillation and classifying sleep stages. Future improvements could involve optimizing feature selection and exploring other machine learning models to further improve classification performance.

---


