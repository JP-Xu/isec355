import os
import numpy as np
from czifile import CziFile
import matplotlib.pyplot as plt
from skimage import filters
from scipy.ndimage import gaussian_filter
import csv

#
# This script can calculated the fluorescence line width by drawing three random vertical lines in a 2D confocal image.
# ostu algorithm is used for mask, and gaussian blur radius (standard deviation) is 2.
#
# Usage:
#       Run this script, all channels (dimensions should be changed based on real images), and results will save to a csv file.


def get_czi_files():
    """ This function returns all czi files under current folder as a list.
    """
    results = []
    files = os.scandir()
    for file in files:
        if file.is_file and file.name.endswith('czi'):
            results.append(file.name)
            
    return results

czi_files = get_czi_files() # czi files under current path.
k = 3  # Number of lines
simga_gauss = 2 # Standard deviation for Gaussian kernel
#channel_index = 0 # Define the channel index will be extracted.
folder_name = os.getcwd().split('/')[-1]
n = 1
results = [] # Result list


print('-'*20)
print(" Params used:\n")
print(" - Number of random lines used: {}".format(k))
print(" - Gaussian blur: {}".format(simga_gauss))
print(" - All channel analysis mode")
#print(" - Channel of czi image: {}\n".format(channel_index))
print('-'*20)
print(" Image processing in progress ...\n")
# Looping through all czi files.
for czi_file in czi_files:
    print(" - Processing image {} ... ".format(czi_file), end="")
    # Open czi_file and save all channel info into image_arrays array.
    with CziFile(czi_file) as czi:
        image_arrays = czi.asarray()
        
    for channel_index in range(image_arrays.shape[1]):
    # Save defined channel to im and squeeze to a 2D array.
        im = image_arrays[0][channel_index].squeeze(axis=2)
        
        # Apply gaussian blur to image. 
        im_gauss = im_clean = gaussian_filter(im, sigma = 2) 
        otsu_threshold = filters.threshold_otsu(im_gauss)
        rand_index = np.random.randint(im_clean.shape[1], size = k)
        
        # Average width
        rand_widths = np.sum(im_clean[:, rand_index] > otsu_threshold, axis=0) * 9.17 / 130 # in microns
        width_mean = np.average(rand_widths)
        with_std = np.std(rand_widths)
        results.append({"id" : n, "file_name" : czi_file, "mean" : width_mean, "std" : with_std, "otsu" : otsu_threshold, "channel": channel_index + 1})
        n += 1
    print(" Finished.")

print("\n Analysis finished, results are saved in the current folder to {}_Con_Results.csv.".format(folder_name[:3]))
print('-'*20)
csv_file = "{}_Con_Results.csv".format(folder_name[:3])
csv_columns = ['id', 'file_name', 'mean', 'std', 'otsu', 'channel']
try:
    with open(csv_file, 'w') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=csv_columns)
        writer.writeheader()
        writer.writerows(results)
except:
    print("Cannot saving to csv file")

    
