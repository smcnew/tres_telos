## Lab notebook telomeres

November 2020

Finished re-running most plates 

- Plate 3 gadph (Nov 2) looked bad, water samples contaminated, failed samples, 
low eff. Will substitue results from July 23, first round, looked very good. 
Comparing plates 4 from 28 July and 10 Nov gadph, repeatibility = 94.5 comparing
raw Cqs 
Plate 1 (october 26th vs. july 17th gadph), repeatability = 92.8

- Did not re-run plate 5 gadph in order to conserve sybr. First attempt (31 july)
looks good. (efficiencies 94-108%; all samples passed qc)

still to do: 
-141 telo samples 
- 8 gadph samples

November 30

Cleaned and processed remaining plates. 
-Switched to B standards for 2020-11-25_telo and _gadph, and 2020-11-24_telo. 
-Probably should exclude 2020-11-12_telo from analysis, because low efficiency. 
All samples from that plate were re-run on 11-25. 

BAD PLATES:
2020-11-12_telo: Low efficiency 
2020-11-02_gadph: Contaminated H20, efficiencies super off. 

GOOD PLATES:
2020-10-26_telo/gadph (Plate 1)
2020-10-29_telo/gadph (Plate 2)
2020-11-02_telo / 2020-07-23_gadph (Plate 3) 
2020-11-10_telo/gadph (Plate 4)
2020-11-24_telo / 2020-07-31_gadph (Plate 5)
2020-11-25_telo/gadph (Redos)
2020-12-4_telo (Final redos)

December 9:
Created a script: "qpcr_summary_compilation.R" to clean and compile qpcr results
Returns a dataset with T/S ratios, 1 per bird. 

December 14: 
Created a script: "adults.R" to compile covariates for adults data set. 
Datasets to be compiled: - 
1. ts_results: telomere results
2. extractions: matching extraction (T_) to band number and capture date 
3. physiology data: matching treatment and phys data to blood sample list
4. nest data: matching fledging success to capture data 

