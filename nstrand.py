#!/usr/bin/env python
# -*- coding: utf-8 -*-

import nest
import numpy.random as rnd

## The random generators and functions defined here enable reproducible results using the same seed
## - as long as the number of virtual processes (ie. nbcpu with --platform=Local) is the same!
## Reproducibility with a different number of virtual processes seems impossible with nest (but do we need it?)

def set_seed(nest_seed, python_seed):
  nest.SetKernelStatus({'grng_seed' : nest_seed})
  N_vp = nest.GetKernelStatus(['total_num_virtual_procs'])[0]
  nest.SetKernelStatus({'rng_seeds' : range(nest_seed+1, nest_seed+N_vp+1)})
  global pyRngs
  pyRngs = [rnd.RandomState(s) for s in range(python_seed+N_vp+1, python_seed+2*N_vp+1)]
  global pyMasterRng
  pyMasterRng = rnd.RandomState(python_seed)

