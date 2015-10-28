#!/usr/bin/env python

#
# Example 3:
# Calculate the discrete Frechet distance between two mirrored paths
#

import numpy # array()
from mpl_toolkits.mplot3d import Axes3D # projection='3d'
import matplotlib.pyplot as plt # plot()

from DiscreteFrechetDistance import *


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

# Calculate the discrete Frechet distance
FD = DiscreteFrechetDistance(P, Q)
cm = FD.getCouplingMeasure()

# Get the coupling sequence
cm_seq = FD.getCouplingSequence()

# Print the coupling sequence and the distance between points
# Note: points in the sequence are number from 1 to p for path P, and 1 to q for path Q.
print ''
print 'Coupling sequence: '
print '  Path P    Path Q    Distance '
for i in range (0, len(cm_seq)):
    u = P[(cm_seq[i,0]-1),:] # minus one because P is indexed from zero
    v = Q[(cm_seq[i,1]-1),:] # minus one because Q is indexed from zero
    dist = numpy.sqrt(numpy.sum((u-v)**2))
    print '    %2i        %2i       %f ' % (cm_seq[i,0], cm_seq[i,1], dist)

print 'Frechet distance = %f \n\n' % cm

# Plot the two paths
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot(Q[:,0],Q[:,1],Q[:,2], '.-', color='red', linewidth=2, label='Q')
ax.plot(P[:,0],P[:,1],P[:,2], '.-', color='blue', linewidth=2, label='P')
#title_str = "title"
#ax.set_title(title_str)
ax.set_xlim([min(P[:,0])-1, max(P[:,0])+1])
ax.set_ylim([-5, 5])
ax.set_zlim([-1, 3])
ax.set_xlabel('x', labelpad=20)
ax.set_ylabel('y', labelpad=40)
ax.yaxis.set_rotate_label(False)
ax.set_zlabel('z', labelpad=40)
ax.view_init(azim=-145, elev=45)
legend = ax.legend(numpoints=1, loc='best', bbox_to_anchor=(1, 0.5), fancybox=True, shadow=True)

frechet_dist_str = 'Frechet distance = %f' % cm
#ax.text(0,0,0, frechet_dist_str, color='black', fontsize=12)
plt.figtext(0.13,0.07, frechet_dist_str, color='black', fontsize=12)

# Plot the coupling sequence
for i in range (0, len(cm_seq)):
    line_x_values = numpy.array([ P[(cm_seq[i,0]-1),0], Q[(cm_seq[i,1]-1),0] ])
    line_y_values = numpy.array([ P[(cm_seq[i,0]-1),1], Q[(cm_seq[i,1]-1),1] ])
    line_z_values = numpy.array([ P[(cm_seq[i,0]-1),2], Q[(cm_seq[i,1]-1),2] ])
    ax.plot(line_x_values, line_y_values, line_z_values, color=[0.5,0.5,0.5])

plt.show()
