# masld_biomarker_estimates
Supports description of true fibrosis biomarker performance.  All datasets are provided for the scripts to be run independently. 

# script sequence

## 00_nafld_f_distribution
Provides script to generate fibrosis distribution per ground truth fibrosis stage and fibrosis stage probability for given fibrosis area from original data.  Generates manuscript figure 1.

## 01_generate_patient_pop
Script to generate a population representative of patients entering biomarker studies.  Outputs a csv file for use in future scripts.

## 02_sampling_variability_comparison
Provides a comparison of simulated biopsy sampling variability with that described by Ratziu et al for simultaneous biopsy of persons with NAFLD.

## 03_bx_performance
Estimates performance characteristics of simulated liver biopsy (area under the curve, sensitivity and specificity) for the identification of advanced and significant fibrosis.  Generates manuscript figure 2.

## 04_blood_biomarker_true_performance
Solves equations described by Waikar et al (JASN 2012) for "true" biomarker performance for commonly used blood-based biomarkers where the target condition is advanced fibrosis.  The apparent performance of these biomarkers is described by Mozes et al (Gut 2021).  Generates manuscript table 1.

## 05_imaging_biomarker_true_performance
Parallel to 04_blood_biomarker_true_performance for imaging biomarkers of fibrosis.  The target condition is advanced fibrosis.  The apparent performance of these biomarkers is described by Selveraj et al (J Hepatol 2021).  Generates manuscript table 2 and figure 3.