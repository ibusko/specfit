

Notebooks and supporting code in this branch are modified versions from their 
original counterparts in master. Here, they exercise the line finding and 
fitting functionality available in specutils. 

The main notebook is linefit_demo_specutils.ipynb. It uses Kriss' model definition
file (n5548_models_specfit.py), which is a modified version of the original file
(n5548_models.py). The modifications were necessary so it can be accepted by 
specutils' line_fit function. 

The other notebook is specfit_demo_specutils.ipynb, which is a re-cast of the 
original demo notebook in master. It uses astropy modelling and fitting, but
stores data in specutils data structures. The model file is the original one.


Contents of this package:

proto:   prototype code

         linefit_demo_specutils.ipynb  - specutil's line finding and line fitting
         specfit_demo_specutils.ipynb  - astropy's model fitting using specutils data
                                         structures

data:    data files

         n5548/n5548_mean_g130mb4.asc:  data from Trello board (Kriss NGC 5548 fitting)
         n5548/n5548_lyalpha_sample.dat: regions defined for data above
