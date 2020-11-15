import numpy as np
import matplotlib.pyplot as plt
from time import time
from utility import *
from sim import *

def simulate(params,sim_type="nucleoid"):
    """ Driver to perform the simulation
    - params: dictionary contains parameters for simulation
    - sim_type: string that can be one of the following: "nucleoid" "free","scaling","rg".
        "nucleoid": simulates the polymer in a confined nucleoid space
        "free": simulates the polymer in free unbound space
        "scaling": calculates the scaling between polymer size and number of segments
        "rg": calculates the radius of gyration of polymer simulated
    """
    assert sim_type.lower() == "nucleoid" or sim_type.lower() == "free" or sim_type.lower() == "scaling" or sim_type.lower() == "rg", "Invalid simulation type, use sim_type keywords: nucleoid, free, Rg or scaling"

    if sim_type.lower() == "nucleoid":
        assert "b" in params, "b must present in params"
        assert "n" in params, "n must present in params"
        b = params["b"]
        n = params["n"]
        if "c" in params:
            c = params["c"]
        else:
            c = 7
        if "solvent" in params:
            solvent = params["solvent"]
        else:
            solvent = "ideal"
        if "trial_limit" not in params:
            trial_limit = 1000
        start = time()
        res = sim_bound(b,n,c,solvent,trial_limit)
        print("Simulated %d mg/ml nucleoid in %s solvent, time elapsed %.2f seconds \n" %(c,solvent,time()-start))
        return res

    elif sim_type.lower() == "free":
        assert "b" in params, "b must present in params"
        assert "n" in params, "n must present in params"
        b = params["b"]
        n = params["n"]
        if "solvent" in params:
            solvent = params["solvent"]
        else:
            solvent = "ideal"
        start = time()
        res = sim_free(b,n,solvent)
        print("Simulated a free chromosome in %s solvent, time elapsed %.2f seconds \n" %(solvent,time()-start))
        return res

    elif sim_type.lower() == "scaling":
        if "ns" in params:
            ns = params["ns"]
        else:
            ns = np.unique(np.arange(10,50,2).astype("int"))
        if "angle" in params:
            rangle = params["angle"]
        else:
            rangle = (0,np.pi)
        if "repeat" in params:
            repeat = params["repeat"]
        else:
            repeat = 1000
        return sim_scaling(ns,rangle,repeat)


    elif sim_type.lower() == "rg":
        assert "b" in params, "b must present in params"
        assert "n" in params, "n must present in params"
        b = params["b"]
        n = params["n"]
        if "solvent" in params:
            solvent = params["solvent"]
        else:
            solvent = "ideal"
        if "repeat" in params:
            repeat = params["repeat"]
        else:
            repeat = 100

        start = time()
        res = sim_rg(b,n,solvent,repeat)
        print("Calculated %d samples of radius of gyration of chromosomes in %s solvent, time elased %.2f seconds \n" %(repeat,solvent,time()-start))
        return res
