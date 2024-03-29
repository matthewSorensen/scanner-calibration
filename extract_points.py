import click
import numpy as np
import sys

from skimage import io, measure, morphology
from skimage.color import rgb2gray
from skimage.filters import threshold_otsu, gaussian
from skimage.segmentation import flood

import matplotlib.pyplot as plt

def contained(inner,outer):
    return outer[0] < inner[0] and outer[1] < inner[1] and inner[2] < outer[2] and inner[3] < outer[3]

def identify_crop(binary_image, raw_image, min_size = 100): 
    labels = measure.label(morphology.binary_dilation(binary_image), background = True)
    
    good_regions = []
    for p in measure.regionprops(labels):
        xmin,ymin,xmax,ymax = p.bbox
        if  (xmax - xmin) < min_size or (ymax - ymin) < min_size:
            continue
        good_regions.append(p)
              
    if len(good_regions) == 0:
        print("Unable to find any crop region - check thresholding parameters and minimum size and try again")
        return None
    elif len(good_regions) > 1:
        print("Found too many crop regions - check thresholding parameters, decrease minimum size, or manually crop image.")
        return None

    # Now we mask out the raw image
    p = good_regions[0]
    xmin,ymin,xmax,ymax = p.bbox

    m = 5
    mask = morphology.binary_erosion(p.image_filled, footprint=[(np.ones((m, 1)), 1), (np.ones((1, m)), 1)])
    crop = raw_image[xmin:xmax, ymin:ymax]

    crop[np.logical_not(mask)] = sum(sum(crop * mask)) / sum(sum(mask))
    return np.array([xmin, ymin]), crop
    
def identify_roi(binary_image, n):
    
    labels = measure.label(morphology.binary_erosion(binary_image))
    m, expected= labels.max(), n*n + 3
    if m < expected:
        print("Too few regions of interest found - check thresholding parameters with --plot")
        return None
    elif 2 * expected < m:
        print("Too many regions of interest found - check thresholding parameters with --plot")
        return None
    
    boxes = np.zeros((m,4),dtype = int)
    boxes[:,0:2] = binary_image.shape
    
    a,b = np.nonzero(labels)
    for i,x in enumerate(a):
        y = b[i]
        l = labels[x,y] - 1
        bbox = boxes[l]
        boxes[l,:] = min(x,bbox[0]),min(y,bbox[1]),max(x,bbox[2]),max(y,bbox[3])
        
    sizes = (boxes[:,2] - boxes[:,0]) * (boxes[:,3] - boxes[:,1])
    idx = np.argsort(0 - sizes)
    circles, marks = list(idx[0:3]), idx[3:3 + n*n]
    
    circled = np.zeros_like(marks)
    for i,j in enumerate(marks):
        if not circles:
            break
        for k,c in enumerate(circles):
                if contained(boxes[j],boxes[c]):
                    circled[i] = 1
                    circles.pop(k)
                    break
    
    return labels, boxes[marks], marks + 1, circled

def find_mark_location(image, binary, box):
    a,b,c,d = box
    offset = int(0.25 * (c - a))
    roi = image[a-offset:c+offset,b-offset:d+offset].astype(float)
    binary = binary[a-offset:c+offset,b-offset:d+offset]
    center = np.unravel_index(np.argmax(roi), roi.shape)

    masked = flood(binary, center) * roi
    total = np.sum(masked)

    x = np.sum(np.arange(0,masked.shape[0]) @ masked) / total
    y = np.sum(masked @ np.arange(0,masked.shape[1])) / total

    return  np.array([x,y]) + np.array([a - offset, b - offset])

@click.command()
@click.option('--negative/--positive', default = True, help = "Default case (--negative) is light marks on a dark background. --positive inverts this.")
@click.option('--crop/--no-crop', default = True, help = "Attempt to crop a full scan to the boundaries of a test coupon")
@click.option('--threshold', default = 2.5, type = click.FloatRange(0.0, 20), help = "User-adjustable parameter for image binarization - larger values bias towards black.")
@click.option('--points', default = 11, type = int, help = "Number of points on the side of the test pattern")
@click.option('--plot/--silent', default = False, help = "Display the segmented image and final result, or just output points.")
@click.argument('input', type=click.File('rb'), nargs = -1)

def cal(*args, **kwargs):
    if kwargs['crop'] and not kwargs['negative']:
        print("Cropping is not supported in negated mode")
        return
    
    print("X, Y, Fiducial, File")
    for x in kwargs['input']:
        run_extraction(x,args, kwargs)
    
def run_extraction(fp,args,kwargs):
    
    # Load the image, convert it to greyscale, and then binarize it
    color_image = io.imread(fp)
    image = rgb2gray(color_image)
    smoothed = gaussian(image, sigma=2)
    thresh = kwargs['threshold'] * threshold_otsu(smoothed)
    binary = (smoothed > thresh)

    offset = np.zeros(2)
    if kwargs['crop']:
        crop = identify_crop(binary, image)
        if crop is None:
            return
        offset, image = crop
        smoothed = gaussian(image, sigma=2)
        thresh = kwargs['threshold'] * threshold_otsu(smoothed)
        binary = (smoothed > thresh)
            
    
    if not kwargs['negative']:
        binary = np.logical_not(binary)
        
    if kwargs['plot']:
        plt.imshow(binary)
        plt.show()

    # Find the regions containing the marks
    roi = identify_roi(binary, kwargs['points'])
    if roi is None:
        return
    
    labels, bboxes, mark_labels, circled = roi
    coords = np.array(list(find_mark_location(smoothed,binary, box) for box in bboxes))
        
    if kwargs['plot']:
        plt.imshow(smoothed)
        plt.scatter(coords[:,1], coords[:,0], c = circled, cmap = 'GnBu')
        plt.show()

    # Then print the coordinates to stdio as a CSV
    for i in np.argsort(0 - circled):
        y,x = offset + coords[i]
        print(f"{x}, {y}, {circled[i] == 1}, {fp.name}")

if __name__ == '__main__':
    cal()

