##  Project: DB/DB
### Study ID: 
### Scientist: David Cheng
### Data Analysis: Davit Sargsyan 
### Created: 09/19/2017 

---    

## Table of Contents
[File Legend](#files)  
[Daily Logs](#logs)   
[Results](#results)   

## File Legend<a name="files"></a>
### Scripts
**/source/dbdb_methylseq_analysis_v3.R**: current (05/15/2018) version of the script for Kidney RO1 Figure 7A, B & C    
**/source/dbdb_rnaseq_deseq2_v1.R**: current (05/15/2018) version of RNA-seq data processing and analysis with DESeq2    
**/source/dbdb_rnaseq_analysis_v2.R**: current (05/15/2018) version of the script for Kidney RO1 Figure 7 D. Uses data created by **dbdb_rnaseq_deseq2_v1.R** script    

### Processed Data
**/data/methyl_seq/combined_DD.csv**: current (05/15/2018) methyl-seq count data; used for Kidney RO1 Figure 7A, B & C 
**/data/methyl_seq/combined.dedup.csv**: data from old samples. No longer used as it has very low horizontal coverage    
**/data/rna_seq/combraw.count**: current (05/15/2018) RNA-seq expression data; used for Kidney RO1 Figure 7D     

### Sample Legends
**/docs/dbdb_David_rerun_dec2017.xlsx**: current (05/15/2018) methyl-seq sample legend.    
**/doc/dbdb_methylseq_legend.xlsx**: FastQ file legend. Where is the legend for RNA-seq files?    
    
### Original FastQ files are located here (4Tb internal drive):    
**/media/administrator/datastorage/2017/Methyl_seq/July/pools/Kong_MouseMethyl_pool5.zip**       
**/media/administrator/datastorage/2017/Methyl_seq/July/pools/Kong_MouseMethyl_pool6.zip**     
   
### All alignemnt files and documents are located here (4Tb internal drive):    
*/media/administrator/datastorage/Processed_BAM_Files/David_MethylSeq_Processed/BAM_Files*   
NOTE: Redo realignment, get BAM files and deduplicate!   
NOTE: some of the files were too large to push to GitHub; saved on the Lab230 machine in *data* folders.

## Daily Logs<a name="logs"></a>
### 05/18/2018
* New DNA data (combined_DD.csv); using nn10 annotation

### 05/15/2018
* Cleaned the script for Kidney RO1 Figure 7 A, B & C (DNA), and 7 D (RNA)

### 12/06/2017
* Reorganized files: all DB/DB files are moved to a new GitHub repository *dbdb*, except results directly related to Kindey RO1 (kept in the original *kidney.ro1* repository).

### 09/27/2017
* Added mapping to oxidative stress and inflamation pathways

### 09/19/2017
* Started

## Results<a name="results"></a>