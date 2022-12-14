# SARS-CoV-2 RBD web application
Observe the change of RBD affinity to ACE2

#### The raw data of RBD affinity score is from this amazing paper:
Starr, T. N., Greaney, A. J., Hilton, S. K., Ellis, D., Crawford, K. H., Dingens, A. S., ... & Bloom, J. D. (2020). Deep mutational scanning of SARS-CoV-2 receptor binding domain reveals constraints on folding and ACE2 binding. cell, 182(5), 1295-1310.

#### This web application was developed by Python and Flask

----------------

### Overview
1. Upload single or multiple fasta files of SARS-Cov-2 strains 
2. Press the button of "show results"

*The web will automatically calculate the sum of RBD score and output the summary table and plot in a html format for each fasta file in the output folder.*

### Installation
#### Clone source code
```
git clone https://github.com/Lilian0296/SARS-Cov-2-RBD-web-application.git

cd RBD_app

python RBD_score.py

```
### Open web
1. Change the port to 8080
2. http://localhost:8080/

### Notice
1. Fasta files for the test were downloaded from NCBI. (The web only takes in-frame fasta files)
2. If there are mutations in RBD domains which are used for grabbing, the web may not work. You may use this tool to get RBD and RBM information:
<https://github.com/Lilian0296/SARS-Cov-2-Spikeinfo.git>

### References
----------

[1] Starr, T. N., Greaney, A. J., Hilton, S. K., Ellis, D., Crawford, K. H., Dingens, A. S., ... & Bloom, J. D. (2020). Deep mutational scanning of SARS-CoV-2 receptor binding domain reveals constraints on folding and ACE2 binding. cell, 182(5), 1295-1310
<https://doi.org/10.1016/j.cell.2020.08.012>



