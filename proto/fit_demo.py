import sys

import fit_functions


data = {'n5548':  ["../data/n5548/",  "n5548_mean_g130mb4.asc",  "n5548_lyalpha_sample.dat", "sfn5548_lyalpha2"],
        'obs1g_1':["../data/simple/", "obs1g.dat",  None, "obs1g_model_1.dat"],
        'obs1g_2':["../data/simple/", "obs1g.dat",  None, "obs1g_model_2.dat"],
        'obs1gn': ["../data/simple/", "obs1gn.dat", None, "obs1g_model_2.dat"]

}

if __name__ == "__main__":

    fit_functions.process_data(data[sys.argv[1]])

