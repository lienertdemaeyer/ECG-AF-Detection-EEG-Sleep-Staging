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

> **Figure 1:** Example of ECG preprocessing steps showing removal of artifacts.

![Preprocessing](https://github.com/user-attachments/assets/9e7a50ef-7a5a-4d3b-b239-f517b677d534)

---

### Artifact Removal
Various filters were applied to remove the following artifacts:
- **Powerline Interference**: Removed using a notch filter.
- **Motion Artifacts**: Filtered using derivative-based and moving average filters.
- **Baseline Wander**: Handled by applying a low-pass filter.

> **Figure 2:** Filtered ECG signal after artifact removal.

![Artifact Removal](https://github.com/user-attachments/assets/731163b6-e836-4829-8f74-4fceeda77979)

---

### Feature Extraction
We extracted features from the ECG signal, including:
- **Morphological features**: Signal energy, wavelet subband energies.
- **Heart-rate variability (HRV) features**: RMSSD, pRR50, and minRR.

> **Figure 3:** R-peak detection for normal and AF segments. Tachogram shown for both cases.

![R-Peak Detection](https://github.com/user-attachments/assets/7e9ef37e-f921-4dc7-835e-5ad5225dcf2d)

---

### Unsupervised Classification
Using Principal Component Analysis (PCA) and K-means clustering, the ECG segments were classified into 'normal' and 'AF' categories.

---

### Evaluation Metrics
The performance of the classification was evaluated based on:
- **Accuracy**: 70.89%
- **Sensitivity**: 71.57%
- **Specificity**: 62.50%

---

## EEG-based Sleep Staging

### Segmentation
The EEG signal was segmented into 30-second intervals.


### Topographic Maps of EEG Sources

Topographic maps were utilized to visually represent the EEG sources identified during Independent Component Analysis (ICA). These maps provide insights into the spatial distribution of electrical activity across the scalp, allowing us to localize brain regions involved in different stages of sleep or related to specific events such as **K-complexes**.

#### K-complex Source Identification

Through the ICA decomposition, each independent component (IC) was analyzed by overlaying its corresponding topographic map with the EEG signal data and convolving it with a K-complex template. Our analysis identified the **7th ICA component** as the likely source of K-complexes. This was supported by the following observations:
- The 7th component exhibited strong activation in the **frontal and central regions**, particularly across the FP1-M2, FP2-M1, F3-M2, F4-M2, and C3-M2 electrodes.
- These regions are well-documented as being active during K-complexes, which are a hallmark of sleep stage 2 (S2).
- **Matched filtering** confirmed this finding, as the 7th component yielded the highest count of detected K-complexes when convolved with the template signal.

> **Figure 4:** Topographic maps of the EEG sources identified via ICA, with components corresponding to various brain and eye activity regions.

![EEG Segmentation](https://github.com/user-attachments/assets/713cf17d-199a-4ba7-a30a-b11776b52266)

#### The Role of EOG Channels in ICA

In EEG analysis, **EOG (Electrooculogram) channels** can sometimes be included in ICA to account for eye movement artifacts. However, in this study, it was determined that excluding EOG channels yielded cleaner brain activity components, particularly in relation to K-complex detection. The **4th and 5th ICA components** exhibited significant activity linked to eye movements, as indicated by their spatial distribution over the **EOG left and EOG right** electrodes. These channels primarily capture:
- **Eye movement artifacts**, such as blinks or horizontal saccades, which can introduce noise into the decomposition process and obscure brain signals.
  
By excluding these EOG channels, the algorithm was able to focus more effectively on EEG signals related to brain activity, improving the detection clarity and accuracy of sleep-specific events like K-complexes.


---

### K-complex Detection
K-complexes were detected using Independent Component Analysis (ICA) and template matching techniques.

> **Figure 5:** Template matching for K-complex detection using matched filters and cross-correlation.

![K-Complex Detection](https://github.com/user-attachments/assets/486d29e0-829c-4cf6-b466-83f159eb6892)

---

### Power Spectral Density Analysis
The Power Spectral Density (PSD) of the EEG signal was computed for 10 different frequency bands to analyze the signal's frequency characteristics during different sleep stages.

---

### Classification Models
A neural network classifier with two hidden layers was implemented to classify the EEG signal into sleep stages (Wake, S2, S3). Feature selection was performed using both manual methods and Linear Discriminant Analysis (LDA).

---

## Results

### Confusion Matrix for Full Model

> **Figure 6:** Confusion matrix showing the performance of the full model for EEG sleep staging classification.

![Confusion Matrix](https://github.com/user-attachments/assets/f96bb8c5-8c01-41e0-8b5c-0f082d7b964d)

- **ECG Classification Results**: The unsupervised model achieved an accuracy of 70.89% for AF detection.
- **EEG Classification Results**: The neural network classifier achieved an accuracy of 72.79% with the full feature set and 81.80% with a reduced set of features.

---

## Conclusion

This project demonstrates the effectiveness of signal processing techniques combined with machine learning models for biomedical data analysis, specifically for detecting atrial fibrillation and classifying sleep stages. Future improvements could involve optimizing feature selection and exploring other machine learning models to further improve classification performance.

---


