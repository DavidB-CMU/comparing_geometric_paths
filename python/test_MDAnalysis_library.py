#!/usr/bin/env python

#
# Compare the performance of two implementations
# for calculating the discrete Frechet distance:
#
#   the DiscreteFrechetDistance class
#     &
#   the MDAnalysis library
#
# You will need to install the python MDAnalysis library:
# sudo pip install --upgrade MDAnalysis
#


import numpy # array()
import time

from DiscreteFrechetDistance import *

import MDAnalysis.analysis.psa


# Define two paths:

P = numpy.array([[0.0, 1.0, 1.0],
                 [1.0, 1.0, 1.0],
                 [2.0, 1.0, 1.0],
                 [3.0, 1.0, 1.0],
                 [4.0, 1.0, 1.0],
                 [5.0, 0.0, 1.0],
                 [6.0, 1.0, 1.0],
                 [7.0, 1.0, 1.0],
                 [8.0, 1.0, 1.0],
                 [9.0, 1.0, 1.0]]);

Q = numpy.array([[0.0, 2.0, 1.0],
                 [1.0, 2.0, 1.0],
                 [2.0, 2.0, 1.0],
                 [3.0, 2.0, 1.0],
                 [4.0, 2.0, 1.0],
                 [5.0, 3.0, 1.0],
                 [6.0, 2.0, 1.0],
                 [7.0, 2.0, 1.0],
                 [8.0, 2.0, 1.0],
                 [9.0, 2.0, 1.0]]);


# Calculate the discrete Frechet distance:

# Using DiscreteFrechetDistance class:
t0 = time.time()
FD = DiscreteFrechetDistance(P, Q)
cm1 = FD.getCouplingMeasure()
t = (time.time() - t0)
print ''
print 'Using own function'
print 'Frechet distance = ', cm1
print 'elapsed time = ', t
print '\n'

# Using function from the MDAnalysis library
t0 = time.time()
cm2 = MDAnalysis.analysis.psa.discrete_frechet(P, Q)
t = (time.time() - t0)
print 'Using MDAnalysis library'
print 'Frechet distance = ', cm2
print 'elapsed time = ', t
print '\n'
