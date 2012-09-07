#!/usr/bin/env python
# -*- coding: utf-8 -*-

# calibration_TIRF.py
# v0.2dev
# last edited 2012-09-06 by MB

"""
determines local maxima in calibration-slides (fluorescent beads),
retuns txt-file with x- and y-coordinates of 'donor'- and 'acceptor' beads
and their relative distances for IGOR-import

uses Numpy, Scipy, Matplotlib and TIFFfile
"""

import os, datetime
import numpy as np 									
import scipy.ndimage as ndimage
import scipy.ndimage.filters as filters
import TIFFfile
import matplotlib.pyplot as plt

def importTIF(directory):
    """ Import tif-files from directory """
 
    filesTIF = []
    for files in os.listdir(directory):
        if files.endswith('.tif'):
            filesTIF.append(files)
    return filesTIF

def localMaxima(sliceString, threshold_Donor, threshold_Acceptor,
                neighborhood):
    """ Input Tiff-Image as string, Return Array with Maxima-coordinates """
    
    data = TIFFfile.imread(sliceString)
    
    data_Donor, data_Acceptor = data[:,:128], data[:,128:]
    
    
    data_max_Donor = filters.maximum_filter(data_Donor, neighborhood)
    maxima_Donor = (data_Donor == data_max_Donor)
    data_min_Donor = filters.minimum_filter(data_Donor, neighborhood)
    diff_Donor = ((data_max_Donor - data_min_Donor) > threshold_Donor)
    maxima_Donor[diff_Donor == 0] = 0

    labeled_Donor, num_objects_Donor = ndimage.label(maxima_Donor)
    slices_Donor = ndimage.find_objects(labeled_Donor)
    xy_Donor = []
    for dy,dx in slices_Donor:
        x_center, y_center = (dx.start + dx.stop - 1)/2, (dy.start + dy.stop -
                                                                1)/2 
        xy_Donor.append([x_center, y_center])
        
    data_max_Acceptor = filters.maximum_filter(data_Acceptor, neighborhood)
    maxima_Acceptor = (data_Acceptor == data_max_Acceptor)
    data_min_Acceptor = filters.minimum_filter(data_Acceptor, neighborhood)
    diff_Acceptor = ((data_max_Acceptor - data_min_Acceptor) > 
                                                threshold_Acceptor)
    maxima_Acceptor[diff_Acceptor == 0] = 0

    labeled_Acceptor, num_objects_Acceptor = ndimage.label(maxima_Acceptor)
    slices_Acceptor = ndimage.find_objects(labeled_Acceptor)
    xy_Acceptor = []
    for dy,dx in slices_Acceptor:
        x_center, y_center = (dx.start + dx.stop - 1)/2, (dy.start + dy.stop -
                                                                1)/2 
        xy_Acceptor.append([x_center, y_center])
        
    return {"donor": np.array(xy_Donor), "acceptor": np.array(xy_Acceptor)}

def plotMaxima(maximaDict, calibData, correspondingMaxima, cleaned):
    """
    Input dictionary with donor- and acceptor-maxima.
    Also input calibration-image as numpy-array.
    Plot the Donor- and Acceptormaxima.
    """
    
    plt.figure(tifName)
    plt.subplot(121)
    plt.plot(maximaDict['donor'][:,0], maximaDict['donor'][:,1],'go')
    plt.plot(maximaDict['acceptor'][:,0], maximaDict['acceptor'][:,1],'ro')
    plt.plot(correspondingMaxima[:,0], correspondingMaxima[:,1], 'bx', 
             markersize=12)  
    plt.plot(correspondingMaxima[:,2], correspondingMaxima[:,3], 'cx', 
             markersize=12)
    plt.imshow(calibData[:, :128] + calibData[:, 128:], 
               interpolation='nearest')
    plt.subplot(122).get_yaxis().set_visible(False)
    plt.plot(maximaDict['donor'][:,0], maximaDict['donor'][:,1],'go')
    plt.plot(maximaDict['acceptor'][:,0], maximaDict['acceptor'][:,1],'ro')
    plt.plot(cleaned[:,0], cleaned[:,1], 'bx', markersize=12)  
    plt.plot(cleaned[:,2], cleaned[:,3], 'cx', markersize=12)
    plt.imshow(calibData[:, :128] + calibData[:, 128:], 
               interpolation='nearest')
    plt.savefig(tifName.split('.')[0] + ".png")
    plt.show()
    
    return

def createCorrespondList(maximaDict, pixels):
	"""
	Take dictionary donor/acceptor and find maxima that correspond
	to each other. Maxima correspond if they are in each others neighbourhood,
	defined by 'pixels'. Neighbourhood equals a circle with radius 'pixels'.
	"""

	correspondCoords = []
	i = 1
	import numpy as nump	
	for acceptor_xy in maximaDict['acceptor']:
		#print 'Acceptor'
		for donor_xy in maximaDict['donor']:
			#print 'Donor'
			if nump.sqrt(nump.square(acceptor_xy[0] - donor_xy[0]) + 
                           nump.square(acceptor_xy[1] - donor_xy[1])) < pixels:
				print 'Treffer ', i
				i += 1
				correspondCoords.append([donor_xy[0], donor_xy[1], 
                             acceptor_xy[0], acceptor_xy[1], acceptor_xy[0] - 
                             donor_xy[0], acceptor_xy[1] - donor_xy[1]])

	return nump.array(correspondCoords)
	
def cleanUp(orangeRedMaximaArray):
    """
    input array (x_orange, y_orange, x_red, y_red, difx, dify) and
    return same array with Maxima existing more than once removed!
    (will remove every red points being coupled to more than one orange point)
    """
    
    d = {}
    for a in orangeRedMaximaArray:
        d.setdefault(tuple(a[2:4]), []).append(a)
    
    cleaned_red = np.array([v for v in d.itervalues() if len(v) == 1])
    cleaned_red = cleaned_red.reshape(cleaned_red.shape[0], 6)
    
    d = {}
    for a in cleaned_red:
        d.setdefault(tuple(a[0:2]), []).append(a)
        
    cleaned = np.array([v for v in d.itervalues() if len(v) == 1])
    
    return cleaned.reshape(cleaned.shape[0], 6)


# main-program using the declared functions:
if __name__ == '__main__':	
    
    
    now = datetime.datetime.now()
    
    np.savetxt('calib%s_result.txt' %now.strftime("%Y-%m-%d"), 
             np.array(['posx_o\tposy_o\tposx\tposy\tdifx\tdify']), fmt='%s')
             
    directory = raw_input('Directory? Press ENTER for CWD  :')
    if directory == '': directory = '.'
    files = importTIF(directory)    
    
    gesamtMaxima = np.empty((1,6))    
    
    for tifName in files:    
    
     
        calibData = TIFFfile.imread(tifName)
     
        maximaDict = localMaxima(tifName, 300, 300, 5)
    
        correspondingMaxima = createCorrespondList(maximaDict, 10)
        
        cleaned = cleanUp(correspondingMaxima)        
        
        gesamtMaxima = np.append(gesamtMaxima, cleaned, axis=0)        
        
        dataFile = open('calib%s_result.txt' %now.strftime("%Y-%m-%d"), 'a')
        
        np.savetxt(dataFile, cleaned, fmt='%d', delimiter='\t')

        dataFile.close()
        
        plot = raw_input('Wanna plot? [y]/n  :')
    
        if plot == "yes" or plot == "Yes" or plot == "YES" or plot == "y" \
                                                             or plot == "Y":
            plotMaxima(maximaDict, calibData, correspondingMaxima, cleaned)
    
    plotGesamt = raw_input('Gesamtplot? [y]/n :')
    
    gesamtMaxima = np.delete(gesamtMaxima, 0, 0)
    
    if plotGesamt == "yes" or plotGesamt == "Yes" or plotGesamt == "YES" or \
                                plotGesamt == "y" or plotGesamt == "Y":
        plt.figure('Alle gefundenen Maxima')
        plt.plot(gesamtMaxima[:,0], gesamtMaxima[:,1], 'gx', markersize=12)  
        plt.plot(gesamtMaxima[:,2], gesamtMaxima[:,3], 'rx', markersize=12)
        plt.plot(gesamtMaxima[:,0], gesamtMaxima[:,1], 'go')  
        plt.plot(gesamtMaxima[:,2], gesamtMaxima[:,3], 'ro')
        plt.axes().set_aspect('equal')
        plt.savefig('gesamtMaxima.png')        
        plt.show()