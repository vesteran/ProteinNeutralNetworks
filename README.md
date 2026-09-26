## This is the Supplementory Materials Repository for the Final Project for Practical Course in Modelling and Systems Biology

Contributors:
Andrew Vester
Varsha Karikalan
Jerneja Volmajer

The Main Folder contains the following files:
- 'README.md' - This file
- 'Analysis.ipynb' - A Jupyter Notebook including the generation of the figures used in the main paper
- 'NeutralNetworkTestingandReducingNotebook.ipynb' - A notebook used for various tests and to use a replacement table for the first attempt at reducing the alphabet of the structures
- 'NeutranNetworkProtein.py' - The main script of this project. Includes the logic and implementation of the random walks as well as the more algorithmic attempt at reducing the alphabet of the structures. 
  - this file takes 2 inputs. The first is the path to the original pdb file for analysis the second being the path where you want to put results. 


The results for the analysis done on each of the proteins are found in their respective folders. The folders contain the following files:
- '#.pdb' - PDB Files for each accepted step in the random walk for each protein
- 'results.csv' - CSV file containing the results of the random walk including a the location of each accepted random walk step, the qmean score, and the number of random walk attempts before finding that successful sequence
- '<Protein Name>_Reduced' - a folder containing the pdb files for each protein after they have been transformed into only having certain amino acids