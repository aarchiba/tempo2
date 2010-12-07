import sys, os
import subprocess

import numpy as np

def run_tests(datadir, base):
    basefn = os.path.join(datadir,base)
    fitsname = basefn+".fits"
    parname = basefn+".par"
    outfitsname = basefn+".out.fits"
    textname = basefn+".out.text"
    goodname = basefn+".good.text"

    print "Testing %s" % base
    subprocess.check_call(["tempo2",
                           "-gr", "photons",
                           "-f", parname,
                           "-in", fitsname,
                           "-textoutput", textname,
                           "-fitsoutput", outfitsname])

    o = np.loadtxt(textname)
    g = np.loadtxt(goodname)

    mjdd = (o[:,1]-g[:,1])*86400
    pd = ((o[:,2]-g[:,2])+0.5)%1-0.5
    pd = (pd-np.mean(pd)+0.5)%1+np.mean(pd)-0.5

    print "Time difference (s) mean: %g std: %g" % (np.mean(mjdd), np.std(mjdd))
    print "Phase difference mean: %g std: %g" % (np.mean(pd), np.std(pd))
    print

if __name__=='__main__':
    datadir = os.path.dirname(sys.argv[0])

    for f in ["fermi-fake", "chandra-1821"]:
        run_tests(datadir, f)
