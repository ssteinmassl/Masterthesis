# Scientific libraries
import numpy as np

# Graphic libraries

from threeML import *
#from astromodels.xspec.xspec_settings import *
#from astromodels.xspec.factory import *

#useful libraries
import sys
from glob import glob
import copy
import collections
import warnings
warnings.simplefilter('ignore')
#import vapeplot
import pymultinest
#from powerlaw_synchrotron import Synchrotron_PowerLaw
from threeML.io.file_utils import file_existing_and_readable

#import own functions
sys.path.append('/home/simonste/research/xrb/scripts/modules')
from my_functions import *

from mpi4py import MPI

import time

comm = MPI.COMM_WORLD
rank = comm.Get_rank()


program_name = sys.argv[0]
#arguments = sys.argv[1:]

#sets = int(arguments[0])

#print(sets)

def chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i + n]
        
        
pa = '/home/simonste/research/xrb/data/'
data_path = pa + 'input%s.txt'%'r'
rdata = XYLike.from_text_file("rdata", data_path)
rdata.assign_to_source(source_name = 'Xrb_rband_source')
datar = DataList(rdata)

def ModeltoPointSource(band):
    mod = TemplateModel('Xrb_HSESDsizeincluded_%s'%band,log_interp = False)
    mod.K.fix = True
    mod.scale.fix =True
    
    return PointSource('Xrb_%s_source'%band, 0, 0, spectral_shape = mod)

filter_list = ['rband']
point_source_list = ['ps_r']
PointSdict = {}
for h in range(len(filter_list)):
    PointSdict[point_source_list[h]] = ModeltoPointSource(filter_list[h])
my_modelr = Model(PointSdict['ps_r'])

my_modelr.Xrb_rband_source.spectrum.main.Xrb_HSESDsizeincluded_rband.T_sec.prior = Truncated_gaussian(mu= 4105., sigma = 19, lower_bound = 4000, upper_bound = 4400)
#my_modelg.Xrb_gband_source.spectrum.main.Xrb_HSESDsizeincluded_gband.T_sec.fix = True
my_modelr.Xrb_rband_source.spectrum.main.Xrb_HSESDsizeincluded_rband.Dtemp.prior = Truncated_gaussian(mu= 1300., sigma = 60, lower_bound = 1000, upper_bound = 2400)
#my_modelg.Xrb_gband_source.spectrum.main.Xrb_HSESDsizeincluded_gband.Dtemp.fix = True
my_modelr.Xrb_rband_source.spectrum.main.Xrb_HSESDsizeincluded_rband.Dsize.prior = Truncated_gaussian(mu = 0.810, sigma = 0.050 , lower_bound = 0.725, upper_bound = 0.99)
#my_modelg.Xrb_gband_source.spectrum.main.Xrb_HSESDsizeincluded_gband.Dsize.fix = True
my_modelr.Xrb_rband_source.spectrum.main.Xrb_HSESDsizeincluded_rband.distfac = 0.001
my_modelr.Xrb_rband_source.spectrum.main.Xrb_HSESDsizeincluded_rband.distfac.prior = Truncated_gaussian(mu = -0.005, sigma = 0.007,lower_bound = -1.,upper_bound =1.)
my_modelr.Xrb_rband_source.spectrum.main.Xrb_HSESDsizeincluded_rband.Inkl.prior = Truncated_gaussian(mu = 70.0, sigma =1.5 ,lower_bound = 60, upper_bound=78)
my_modelr.Xrb_rband_source.spectrum.main.Xrb_HSESDsizeincluded_rband.HSWidth.prior = Truncated_gaussian(mu = 11, sigma =6 ,lower_bound =5,upper_bound =45)
my_modelr.Xrb_rband_source.spectrum.main.Xrb_HSESDsizeincluded_rband.HSAwidth.prior = Truncated_gaussian(mu = 0.012, sigma =0.007 ,lower_bound =0.001,upper_bound =0.022)
my_modelr.Xrb_rband_source.spectrum.main.Xrb_HSESDsizeincluded_rband.HSTemp.prior = Truncated_gaussian(mu = 4100, sigma =700 ,lower_bound =3000,upper_bound =9000) 

bsMn = BayesianAnalysis(my_modelr, datar)
bsMn.sample_multinest(n_live_points=400,
                    resume=False,
                    importance_nested_sampling=False,
                    verbose=False,
                    chain_name='chains_eff/HSESDsizeinclonlyr0602/full_band')

time.sleep(10)

if rank == 0:

    bsMn.results.write_to('XrbHSESDsizeinclonlyr0602Multinest.fits')
