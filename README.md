

One of the main purposes of this sprint was to lay down the foundations
of a sort of "language" that would be used to express arbitrarily complex,
multi-component spectral models. The model expressions should be easily
modifiable via simple text editing tools. Also, they should be self-
documenting, and easily executable by plain python code. This would
enable users to record and store arbitrary spectral model specifications
in plain text documents.

The best that could be achieved in that regard was the use of python
itself to write the spectral model expressions. That way, a text file
with model specifications can be easily ingested via a simple 'import'
statement.

An example is stored in file proto/n5548_models.py. An example of
use is documented in the specfit_demo_2.ypnb ipython notebook.

GUI support: TBD


Contents of this package:

proto:   prototype code

         To run:
         - cd to the 'proto' directory
         - run the test fit script with the data set name as a command line argument, as in:

           % python fit_demo.py n5548

           or

           % python fit_demo.py obs1g_1

           Note that astropy cannot handle tied parameters yet. Tied parameters in Jerry's model
           are being handled as free parameters for now.

data:    data files

         n5548:  data from Trello board (Kriss NGC 5548 fitting)
         simple: Jeff's simple 1 and 2 components Gaussian fitting.

         See also the readme file in data/simple/


