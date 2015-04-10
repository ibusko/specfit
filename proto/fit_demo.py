import sys

import fit_functions


data = {'n5548':   ["../data/n5548/",  "n5548_mean_g130mb4.asc",  "n5548_lyalpha_sample.dat", "sfn5548_lyalpha_all_first.dat"],
        'obs1g_1': ["../data/simple/", "obs1g.dat",  None, "obs1g_model_1.dat"],
        'obs1g_2': ["../data/simple/", "obs1g.dat",  None, "obs1g_model_2.dat"],
        'obs1gn_1':["../data/simple/", "obs1gn.dat", None, "obs1g_model_1.dat"],
        'obs1gn_2':["../data/simple/", "obs1gn.dat", None, "obs1g_model_2.dat"],
        'obs2g_1': ["../data/simple/", "obs2g.dat",  None, "obs2g_model_1.dat"],
        'obs2g_2': ["../data/simple/", "obs2g.dat",  None, "obs2g_model_2.dat"],
        'obs2gn_1':["../data/simple/", "obs2gn.dat", None, "obs2g_model_1.dat"],
        'obs2gn_2':["../data/simple/", "obs2gn.dat", None, "obs2g_model_2.dat"],

}

if __name__ == "__main__":

    fit_functions.process_data(data[sys.argv[1]])

