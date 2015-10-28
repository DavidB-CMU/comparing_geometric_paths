#!/usr/bin/env python

# Copyright (c) 2015, Carnegie Mellon University
# All rights reserved.
# Authors: David Butterworth <dbworth@cmu.edu>
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# - Redistributions of source code must retain the above copyright notice, this
#   list of conditions and the following disclaimer.
# - Redistributions in binary form must reproduce the above copyright notice,
#   this list of conditions and the following disclaimer in the documentation
#   and/or other materials provided with the distribution.
# - Neither the name of Carnegie Mellon University nor the names of its
#   contributors may be used to endorse or promote products derived from this
#   software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.


# Calculate the discrete Frechet distance (the coupling measure)
# between two paths P and Q.
#
# Usage:
#   dfd = DiscreteFrechetDistance(P, Q)
#   cm = dfd.getCouplingMeasure()
#   cm_seq = dfd.getCouplingSequence()
#
# Input: Two piece-wise linear paths P and Q
#        Arrays where the rows are the waypoints, and columns are 
#        x,y,z values of each waypoint.
#
# Output: cm The coupling measure
#         cm_seq The coupling sequence
#
#
# David Butterworth, 2015
#
# Algorithm by:
# T. Eiter and H. Mannila, "Computing Discrete Frechet Distance",
# Technical Report CD-TR 94/64, Christian Doppler Laboratory for Expert
# Systems, Vienna University of Technology, 1994
#
# Back-tracking algorithm by:
# Zachary Danziger, 2011
#


import numpy


class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'


class DiscreteFrechetDistance(object):
    """
    A class for calculating the discrete Frechet distance.
    """

    def __init__(self, path1, path2):

        # Check input
        if (type(path1) != numpy.ndarray) or (type(path2) != numpy.ndarray):
            print bcolors.FAIL + "ERROR: path1 and path2 must be numpy arrays" + bcolors.ENDC
            return
        if (path1.ndim != 2) or (path2.ndim != 2):
            print bcolors.FAIL + "ERROR: path1 and path2 must be 2-dimensional numpy arrays" + bcolors.ENDC
            return

        (path1_rows, path1_cols) = path1.shape
        (path2_rows, path2_cols) = path2.shape

        if path1_cols != path2_cols:
            print bcolors.FAIL + "ERROR: path1 and path2 must have the same number of columns" + bcolors.ENDC
            return
        elif (path1_cols < 2) or (path2_cols < 2):
            print bcolors.FAIL + "ERROR: path1 and path2 must have 2 or more columns" + bcolors.ENDC
            return

        self.P = path1
        self.Q = path2

        self.p = path1_rows # Number of points in P
        self.q = path2_rows # Number of points in Q

        # ca: array [1..p, 1..q]
        # (all values must be -1 at start)
        self.ca = -1 * numpy.ones([self.p, self.q]);


    def getCouplingMeasure(self):
        """
        Return the coupling measure (the discrete Frechet distance).
        """

        def d(u, v):
            """
            function d(i,j)
            Calculate the Euclidean distance between two points (u_i, v_i) on the paths P and Q.
            (also called L2 or Pythagorean distance)
            """

            return numpy.sqrt(numpy.sum((u-v)**2))

        def c(i,j):
            """
            function c(i,j)
            Calculate the coupling measure between two paths P and Q.
            by iterating over the waypoints i = 0 to (p-1) on path P, and j = 0 to (q-1) on path Q.

            The original algorithm indexes the first point on each path as one, however python indexes arrays from zero.
            """

            if (self.ca[i,j] > -1):
                return self.ca[i,j]
            elif (i == 0) and (j == 0):
                self.ca[i,j] = d(self.P[0,:],self.Q[0,:])
                return self.ca[i,j]
            elif (i > 0) and (j == 0):
                self.ca[i,j] = max( c(i-1,0), d(self.P[i,:],self.Q[0,:]) )
                return self.ca[i,j]
            elif (i == 0) and (j > 0):
                self.ca[i,j] = max( c(0,j-1), d(self.P[0,:],self.Q[j,:]) )
                return self.ca[i,j]
            elif (i > 0) and (j > 0):
                self.ca[i,j] = max( min([c(i-1,j), c(i-1,j-1), c(i,j-1)]), d(self.P[i,:],self.Q[j,:]) )
                return self.ca[i,j]
            else:
                self.ca[i,j] = numpy.inf

        # Return the coupling measure
        return c((self.p - 1), (self.q - 1))


    def getCouplingSequence(self):
        """
        Return the coupling sequence, which is the sequence of steps
        along each path that result in the coupling measure.
        """

        cm_seq = numpy.zeros([self.p+self.q+1, 2])

        # Create matrix where first row and first column are inf,
        # with matrix ca contained within.
        padded_ca = numpy.inf * numpy.ones([self.p+1, self.q+1])
        padded_ca[1:(self.p+1), 1:(self.q+1)] = self.ca

        Pi = self.p
        Qi = self.q
        count = 0
    
        # Iterate from p down to 2, and q down to 2
        while ((Pi != 1) or (Qi != 1)):

            # Step down the gradient of matrix padded_ca
            min_idx = numpy.argmin([padded_ca[Pi-1,Qi], padded_ca[Pi-1,Qi-1], padded_ca[Pi,Qi-1]])
            if (min_idx == 0):
                cm_seq[count,:] = [Pi-1, Qi]
                Pi = Pi - 1
            elif (min_idx == 1):
                cm_seq[count,:] = [Pi-1, Qi-1]
                Pi = Pi - 1
                Qi = Qi - 1
            elif (min_idx == 2):
                cm_seq[count,:] = [Pi, Qi-1]
                Qi = Qi - 1

            count = count + 1

        # Unlike matlab version, we don't need to subtract 1 (the padding value)
        # from the indices in cm_seq because python indexes arrays from zero.

        # Find the index of the last non-zero value in cm_seq
        last_value_idx = 0
        for i in range (0, len(cm_seq)):
            if (cm_seq[i,0] == 0):
                last_value_idx = i - 1
                break

        # Get the non-zero values
        cm_seq = cm_seq[0:(last_value_idx+1),:]

        # Flip order of rows from bottom-to-top
        cm_seq = cm_seq[::-1]

        # Add the last point of P and Q
        cm_seq = numpy.vstack((cm_seq, [self.p, self.q]))

        return cm_seq
