# ML-to-predict-COI-sub-errors
Towards High-Throughput Biodiversity Analysis: Applying Machine Learning to Predict Whether COI Sequences Contain Substitution Errors


Step 1:
This folder contains the scripts/datasets used to prepare the COI sequences for the ML pipelines. Data preparation for the three animal phyla are divided into three R scripts (“Annelida Data. R”; “Echinodermata Data.R”; “Tardigrada Data.R”). The datasets acquired from BOLD (The Barcode of Life Data Systems) can also be found in this folder (“Annelida Data From BOLD.txt”; “Echinodermata Data From BOLD.txt”; “Tardigrada Data From BOLD.txt”). Within the R scripts, COI sequences get made into FASTA files (on three occasions) to add artificial substitution errors to half the sequences (at 2%, 5% and 10%). These artificial substitution errors are added using three python scripts available in this folder (“Error Generator -2%.py”, “Error Generator -5%.py”, “Error Generator -10%.py”). Once these substitution errors are added, the FASTA files are read back into the R script to continue the data preparation process.

Step 2:
This folder contains the output files from “Step 1.” There are three folders for the different animal phyla. Each folder contains six CSV files containing one-hot encoded COI sequences and their label (1= error; 2= clean). Within each folder three of the files contain the major order sequences at the different errors rates; and the other three files contains the minor order sequences at the different error rates. 

Step 3:
This folder contains the pipelines and scripts used to create the ML models and evaluate their performance. In total there are 54 pipelines to run. Each of the datasets from “Step 2” is run on three different ML pipelines (SVM, LR, RF). When training/validation/testing is done on major orders, only one dataset is used. When training/validation is done on major orders and testing is done on minor orders, two datasets are used. 

Step 4:
This folder contains sample output folders (and their files) after running two pipelines for the phylum Tardigrada at 5% error rate. The first five folders (contains “2ds” in the folder name) are the output folders after using both major and minor orders (2 datasets). The final five folders are the output folders after using just major orders (1 dataset). Each pipeline (from “Step 3”) generates these 5 folders. 
