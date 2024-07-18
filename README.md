#deepTAD: Identifying Topologically Associating Domains  Based on Deep Learning Methods

# About deepTAD

deepTAD, which is based on CNN-transformer model for feature extraction from Hi-C contact matrix, predicting the boundaries in the matrix of new samples by training the model, and subsequently, screening out the false-positive boundaries by using the Wilcox rank-sum test, and finally, completing the TAD assembly by cosine similarity.

# Requirements

* python 3.7.12, tensorflow2.9.0 , numpy 1.21.6, pandas 1.3.5,  scipy 1.7.3, matplotlib  3.7.2.



# Usage <a name="usage"></a>

## First--prepare_data
1. **_You can use "sh download_raw_hic.sh" to download the raw HiC data needed for the experiment.Then, use "sh dump_data_from_hic.sh to generate the contact frequency matrix for different normalization methods (KR or VC),all required files are in the prepare_data folder._**

2. ***Normalized data for the other five cell lines (IMR90, K562, NHEK, HMEC, HUVEC) run the command "sh Cells_observed.sh"***

   


## Second--Sample generation
1. ***Identify the TAD on HIC002 using the six TAD callers 'CaTCH', 'CHAC', 'deDoc', 'DI', 'TopDom', 'Arrowhead', and find the boundaries common to any two methods***

   

2. ***Generate 10×10 sample data from large Hi-C contact matrix matrices based on the generated boundaries and label them to generate training and validation datasets***

   

3. ***The code needed for this section is in the "samples process.ipynb" in the sample generation folder.***

   


## Third --Model

 ***Load the model weights and make predictions to get the TAD boundaries,Model weights for "Model/all.27.h5", run code in "Model/cnn+transformer.ipynb"***



## Final --Filtering false-positive borders and assembling merged TADs

1. ***Filter false positive boundaries with "bash assemble_filterTAD.sh"***

   

2. ***The code to merge TAD is in "filter boundary and merge TAD/mergeTAD.ipynb"***







