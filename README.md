# introgression
The simulation can be run from the file **introgression.f90**. The parameters for the simulation are set in the file **input.in**. The files **seed.in** and **forbiden-sites.in** are also necessary for the simulation. The code is prepared to be run in a loop. All output files will be stored in an individual folder "run_xx" for each realization. The outputs are:

   - pop-new.dat  
   _saves the simulation parameters and genomes (nuclear and mitochondrial) and positions of all individuals of the population at the final generation_
   
   - number0.dat  
   _contains the number of species at each generation_

   - count-mut.dat  
   _saves the number of alleles ''1'' in each segment of each genome every generation; also counts the number of switches 1-0 and 0-1, for the case that the incompatibilities are arranged at random at the start of the simulation; the columns are in this order: sum(mutmit1), sum(mutmit2), sum(mutnuc1), sum(mutnuc2), sum(mutnucneutral), sub10, sub01_

   - survivals.dat  
   _saves the index of the individual from the first generation, from whom the final population inherited the mitochondrial genome_
  
The plots of the manuscript were based on an ensemble of 50 runs for each value of mito-nuclear selection strenght. The data analysis made in Jupyter notebooks (Python) are contained in the files:

   - mutation.ipynb   
   _plots when changing the mitochondrial mutation rate (μ<sub>m</sub>) (fig. 2, 3, and 6a)_

   - selection.ipynb   
   _plots when changing the strength of selection (S) (fig. 4 and 6b)_
   
   - incompatibilities.ipynb   
   _plots when changing the fraction of incompatible interacting genes (f) (fig. 5 and 6c)_
   
  
  REFERENCE: _"Nuclear compensatory evolution driven by mito-nuclear incompatibilities"_ Débora Princepe, Marcus A. M. de Aguiar. arXiv:2403.10411 (2024) https://arxiv.org/abs/2403.10411
