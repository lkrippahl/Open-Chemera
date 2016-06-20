# -*- coding: utf-8 -*-
"""Rotation optimization

This module implements a greedy stochastic hill climbing local search algorithm to 
optimize the sample of orientations used in the rotation search in order to
better adapt it to the irregular shape of the probe protein

"""

import time
import argparse
import numpy as np
import xml.etree.ElementTree as ET
import matplotlib.pyplot as plt

__author__ = "Goncalo Sousa Mendes and Ludwig Krippahl"
__license__ = "Public Domain"
__version__ = "0.1"
__status__ = "Prototype"

NEIGHBOURS = 20#498
"""Size of neighbourhood to check distances"""

SCALE = 1.0
"""Scale of random perturbation for new quaternions"""

ADAPTIVE=1.0
"""Scale change on successful improvement or failure
(1.0 is no change)
"""

MINSCALE=0.01
"""Minumum perturbation scale for adaptive scaling
(never used if ADAPTIVE is 1.0)
"""

RECALC_LIMIT = 1.1
"""Radius variation above which neighbourhood is recomputed"""

KTILES=10
"""Number of percentiles for statistics"""


def read_xml_file(file_name):
    """Read quatrnions from the specific file, 
    Return an array of quaternions, size ((nQuat,4)) 
    The file name must include the path from the current folder
    Param: the file name
    Return: a matrix with the quaternions
    """

    tree = ET.parse(file_name)
    root = tree.getroot()
    quaternions = []
    rotations = root.find('DockRun').find('Rotations')    
    for q in rotations.iter('Rotation'):
        quaternions.append([ q.get('r') , q.get('i')  ,  q.get('j') ,  q.get('k') ])
    return np.array(quaternions,dtype=np.float)
    

def read_points(file_name):
    """Read the points that represent the protein from a file .pts
    The file name must include the path from the python folder
    Param: the file name
    Return: a matrix with the points
    """
    t = list(file_name.split('/'))
    t = t[1].split('-')
    file_name = t[0]
    points = np.zeros((20,3))
    i = 0
    f = open('points/' + file_name + '.pts','r')
    for line in f.xreadlines():
        l = list(line.split(';'))
        for s in range(len(l)):
            points[i,s] = l[s] 
            #print l[s]
        i = i + 1
    f.close()
    return points
    

def write_xml(file_name, new_file_name , quat):
    """Write to a xml file the new sets of quaternions
    The file name must include the path from the python folder
    Param: the file name of the original quaternions
    Param: the file new file name 
    Param: the quaternions to write
    """

    tree = ET.parse(file_name)
    root = tree.getroot()
    rotations = root.find('DockRun').find('Rotations')
    i = 0
    for q in rotations.iter('Rotation'):
        q.set('r', str(quat[i,0]) )
        q.set('i', str(quat[i,1]) )
        q.set('j', str(quat[i,2]) )
        q.set('k', str(quat[i,3]) )
        q.set('ix', str(0))
        i = i+1    
    tree.write(new_file_name)
    

def multiply_quaternions(quats1,quats2):
    """Return a matrix (#quaternions,4) resulting from the pairwise 
    multiplication of quaternions. Matrices quats1 and quats2 must have
    same number of lines and 4 columns
    Param: a set of quaternions
    Param: a set of quaternions
    Return: a matrix with the results
    """
    w1 = quats1[:,0]
    x1 = quats1[:,1]
    y1 = quats1[:,2]
    z1 = quats1[:,3]

    w2 = quats2[:,0]
    x2 = quats2[:,1]
    y2 = quats2[:,2]
    z2 = quats2[:,3]

    res = np.zeros((quats1.shape[0],4))
    
    res[:,0] = w1 * w2 - x1 * x2 - y1 * y2 - z1 * z2
    res[:,1] = w1 * x2 + x1 * w2 + y1 * z2 - z1 * y2
    res[:,2] = w1 * y2 + y1 * w2 + z1 * x2 - x1 * z2
    res[:,3] = w1 * z2 + z1 * w2 + x1 * y2 - y1 * x2
    return res
    
    
def conjugate(quats):
    """Return conjugates of quaternions
    Param: a set of quaternions
    Return: a matrix with the conjugates
    """    
    res = np.zeros(quats.shape)
    res[:,0]=quats[:,0]
    res[:,1]=-quats[:,1]
    res[:,2]=-quats[:,2]
    res[:,3]=-quats[:,3]
    
    return res
    
    
def normalize_quat(quat):   
    """Return the normalized quaternions
    Param: a set of quaternions
    Return: a matrix with the normalized quaternions
    """

    mag = np.sum(quat**2,axis=1, keepdims=True)
    r = quat / np.sqrt(mag)   
    return r


def rotate_points(points,quaternions):
    """Return matrix with the rotated points
    points is matrix of points lines and 3 columns
    quaternions is a matrix with quaternion lines and 4 columns
    result is a matrix with (#quaternions, #points, 3) dimensions
    Param: a set of points
    Param: a set of quaternions
    Return: a matrix with the rotated points
    """    

    res = np.zeros((quaternions.shape[0],points.shape[0],4))  
    res[:,:,1:] = points
    conjugates = conjugate(quaternions) 
    for ix in range(len(points)):
        res[:,ix,:] = multiply_quaternions(quaternions,res[:,ix,:])
        res[:,ix,:] = multiply_quaternions(res[:,ix,:],conjugates)
    return res[:,:,1:]


def point_dists_one_point(base_sets,new_sets):
    """Return vector of max distances between the base_point
    and the new_point (n quaternions to one quaternion)
    Param: a set of base points
    Param: a set of new points
    Return: a matrix with the distances
    """

    diffs = base_sets-new_sets
    dists = np.sum(np.square(diffs),axis=-1)        
    return np.sqrt(np.max(dists,axis=-1))


def quaternion_step(quaternions, neighbours, scale):
    """Create a new set of quaternions with a step of scale relative to
    the difference to the closest quaternion
    Param: a set of quaternions
    Param: a set with the neighbours
    Param: the scales
    Return: a matrix with the new quaternions
    """

    randoms = np.random.uniform(-1,1,(quaternions.shape[0],4))*scale[:,np.newaxis]
    neighs = quaternions[neighbours[:,1],:]
    angle =  1-np.square(np.sum(quaternions*neighs,axis=-1,keepdims=True))
    angle = angle+1.0/quaternions.shape[0]
    new = quaternions + randoms*angle
    return normalize_quat(new)


def compute_neighbours_quats(quaternions,n):
    """Return the indicies of the n closest quaternions from the given quaternion
    Param: a set of quaternions
    Param: the number of neighbours
    Return: a matrix with the indicies
    """

    dists =  1-np.square(np.sum(quaternions*quaternions[n,:],axis=-1))
    ixs = np.argsort(dists)[0:NEIGHBOURS+1]
    ixs = ixs[ixs!=n]
    return ixs

def hill_climbing(quaternions, points, rot_points, number_it):
    """Optimize rotations by hill climbing
    Param: a set of quaternions
    Param: a set of points
    Param: a set with rotate points with th initial quaternions
    Param: the number of iterations
    Return: the decils with the distances to plot the graph
    Return: the averege of updates
    Return: the new set of quaternions
    """    

    it = 0
    # to save the incidies of the closest quaternions to each quaternion
    neighbour_ixs = np.zeros((quaternions.shape[0],NEIGHBOURS),dtype=np.int)
    # to save the rotate
    #rot_points = rotate_points(points, quaternions)
    # to save the decils with the distances to plot the graph
    ktiles = np.zeros((number_it,KTILES+1))
    # to save the minimum distances for each quaternion
    mins = np.zeros((quaternions.shape[0]))
    # to save the max distance in the neighborhood, to control 
    # the need to recalculate the neighborhood
    neighbour_rads = np.zeros((quaternions.shape[0]))
    # to save the scales
    scales = np.ones((quaternions.shape[0]))*SCALE
    # to save the total updates
    total_updates = 0
    
    # initial loop to discover the neighbours of each quaternion
    for n in range (0,quaternions.shape[0]):
        # find the indicies of the closest quaternions
        neighbour_ixs[n,:]=compute_neighbours_quats(quaternions,n)
        # get the rotate points of the closest quaternions
        neighs = rot_points[neighbour_ixs[n,:],:]
        # find the distances of the quaternion to all of its neighbours
        dists = point_dists_one_point(neighs, rot_points[n])
        # update the max distance
        neighbour_rads[n] = np.max(dists)
        # update the min distance
        mins[n] = np.min(dists)    
            
    first_min = np.min(mins)
    print 'first min: ', first_min, np.max(mins)
    start_time = time.time()
    
    
    #main loop that optimizes the quaternions
    while ( it < number_it):
        #  save the number of iterations of the current loop
        updates = 0
        # number of neighborhood recalculations
        recalcs = 0
        # violation of neighborhood
        total_neigh_violation = 0.0
        # to save the minimum distances for each quaternion
        mins = np.zeros((quaternions.shape[0]))
        
        scales[scales<MINSCALE]=MINSCALE
        # create the new set of quaternions
        news_quats = quaternion_step(quaternions, neighbour_ixs, scales)
        # rotate the points
        new_rot_points = rotate_points(points, news_quats)
        # loop that iterates through all the quaternions and test the new one
        for n in range (0,quaternions.shape[0]):
            # get the indicies of the neighbours
            neighs = rot_points[neighbour_ixs[n,:],:]
            # find the distances of the quaternion to all of its neighbours using the old quaternion
            all_dist_old = point_dists_one_point(neighs, rot_points[n])
            old_quat_dist = np.min(all_dist_old)            
            # find the distances of the quaternion to all of its neighbours using the new quaternion
            all_dist_new = point_dists_one_point(neighs,new_rot_points[n])
            new_quat_dist = np.min(all_dist_new)  
            #compare the two distances
            if new_quat_dist>old_quat_dist:
                # the new distance is a good one
                # update the scales
                scales[n]=scales[n]*ADAPTIVE
                # update the rotated points
                rot_points[n] = new_rot_points[n]
                # update the quaternion
                quaternions[n] = news_quats[n]
                # update the mins
                mins[n] = new_quat_dist
                # increase the number of updates
                updates+=1
                # to compare with the old max distance to test the neighborhood
                radius =  np.max(all_dist_new)
            else:
                # the new distance is a bad one (no substitution)
                # update the mins
                mins[n] = old_quat_dist
                # update the scales
                scales[n]=scales[n]/ADAPTIVE
                # to compare with the old max distance to test the neighborhood
                radius =  np.max(all_dist_old)
                
            # test to see if is necessary to recalculate the neighborhood
            if radius > neighbour_rads[n]*RECALC_LIMIT:
                # update the scales
                scales[n]=scales[n]/ADAPTIVE
                # update the total_neigh_violation
                total_neigh_violation+=radius - neighbour_rads[n]
                # find the new neighbours
                neighbour_ixs[n,:]=compute_neighbours_quats(quaternions,n)
                # get its rotated points
                neighs = rot_points[neighbour_ixs[n,:],:]
                # calculate the distances
                dists = point_dists_one_point(neighs, rot_points[n])
                # update the max distance
                neighbour_rads[n] = np.max(dists)  
                # update the mins
                mins[n] = np.min(dists)   
                # update the recalcs
                recalcs+=1
              
        # order the mins
        dist_ixs = np.argsort(mins)
        ile = float(len(mins))/KTILES
        # update data to plot the graph
        for ix in range(KTILES):
            ktiles[it,ix]=mins[dist_ixs[int(ile*ix)]]
        ktiles[it,-1]=mins[dist_ixs[-1]]
                
        print ' '.join(['{0:.2f}'.format(v) for v in ktiles[it,:]])
        elapsed = (time.time() - start_time)/60.0
        print '{0}, elapsed: {1:.2f}m, total: {2:.2f}m'.format(it,elapsed,elapsed*number_it/(it+1.0))
        print 'updates: ',updates, min(scales),max(scales)
        print 'recalcs: ',recalcs, total_neigh_violation/(recalcs+1)
        total_updates+=updates        
        it = it + 1
        
    return ktiles, float(total_updates)/it, quaternions

def write_stats(file_name,stats,total_time):
    """Write statics to file
    Param: the name of the file
    Param: the statistics
    Param: the total time
    """

    ktiles, av_updates = stats    
    ofil = open(file_name+'.stats','w')
    ofil.write('Initial dists:\n '+' '.join(['{0:.2f}'.format(v) for v in ktiles[0,:]]))
    ofil.write('\nFinal dists:\n '+' '.join(['{0:.2f}'.format(v) for v in ktiles[-1,:]]))
    ofil.write('Average of {0:.2f} updates in {1:.2f}m\n'.format(av_updates,total_time))
        
    def line(vals):
        return '\t'.join(['{0:.2f}'.format(v) for v in vals])
    table = '\n'.join([line(ktiles[i,:]) for i in range(ktiles.shape[0])])
    ofil.write('min\tmax\tavg\tstd\n')
    ofil.write(table)
    ofil.close()


def run_faster(file_name, number_it,new_file_name ):
    """Run optimization and write to file    
    """
    
    print 'Starting the ', file_name, ' job...'
    start_time = time.time()
    
    print 'Reading points...'
    points = read_points(file_name)  
    quat = read_xml_file(file_name)
    print quat.shape[0]
    print 'Rotating points...'
    start_time_rotation = time.time()
    rots = rotate_points(points, quat)
    print 'Time for rotation: ', time.time() - start_time_rotation

    stats = hill_climbing(quat, points,rots, number_it)
    new_quat=stats[-1]
    total_time = (time.time()-start_time)/60.0
    write_stats(file_name,stats[:-1],total_time)
    write_xml(file_name, new_file_name,new_quat)
    
   
    print "\nElapsed time:",total_time
    plt.figure()
    plt.plot(range(stats[0].shape[0]),stats[0])
    plt.savefig(file_name+'-'+str(NEIGHBOURS)+'-'+str(SCALE)+'.png',dpi=600)
    plt.close()


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Rotation optimization by hil climbing.')
    
    parser.add_argument('in_file',  action='store', help='Input file, mandatory')
    
    parser.add_argument('--out',  action='store', dest='out_file', default=None,
                        help='Output file (optional)')    
    
    parser.add_argument('--it',  action='store', dest='number_it', default=100,
                        help='number of iterations, default: 100 (optional)', type=int)
    
    
    args = parser.parse_args()
    
    run_faster(args.in_file, args.number_it,
               args.out_file)

