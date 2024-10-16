import numpy as np
from astropy import units as u
from astropy.io import fits
from astropy.nddata import CCDData, Cutout2D
from ccdproc import ImageFileCollection, subtract_bias

import argparse
import logging as log

# Logging
parser = argparse.ArgumentParser(
                    prog='quicktrim',
                    description='quicktrim bias subtracts and trims repeated overscan regions of a ccd.',
                    epilog='This was written to be used with STA1600 CCD')

parser.add_argument('--verbose', '-v', action='store_true')

args = parser.parse_args()
if args.verbose:
    log.basicConfig(format="%(levelname)s: %(message)s", level=log.DEBUG)
    log.info("Verbose output.")
else:
    log.basicConfig(format="%(levelname)s: %(message)s")


# Overscan Calculations
# horizontal
horizontal_taps = 8
tap_width = 1320
overscan_width = 180 # per tap
tap_height = 5280
overscan_height = 20 # per tap
xstart = []
ystart = []
xend = []
yend = []

log.info("Horizontal taps:")
for htap in range(1, horizontal_taps+1):
    log.info(f"{htap}\t {(htap-1)*(tap_width+overscan_width)} \t {(htap-1)*(tap_width+overscan_width)+tap_width}")
    xstart.append((htap-1)*(tap_width+overscan_width))
    xend.append((htap-1)*(tap_width+overscan_width)+tap_width)

# vertical
log.info("Vertical taps:")
vertical_taps = 2
for vtap in range(1, vertical_taps+1):
    log.info(f"{vtap}\t {(vtap-1)*(tap_height+2*overscan_height)} \t {(vtap-1)*(tap_height+2*overscan_height)+tap_height}")
    ystart.append((vtap-1)*(tap_height+2*overscan_height))
    yend.append((vtap-1)*(tap_height+2*overscan_height)+tap_height)

# create slice list to remove overscan regions
# slices will contain data so delete function will use not operator ~ to keep data slices
xdata_slice_list = ()
ydata_slice_list = ()

for (start, end) in zip(xstart, xend):
    xdata_slice_list = xdata_slice_list + (slice(start, end),)

for (start, end) in zip(ystart, yend):
    ydata_slice_list = ydata_slice_list + (slice(start, end),)

log.info('x slices')
for xdata in xdata_slice_list:
    log.info('\t'+str(xdata))
log.info('y slices')
for ydata in ydata_slice_list:
    log.info('\t'+str(ydata))

# get image list
ifc = ImageFileCollection('./data/')

# load bias data
with fits.open('./bias.fits') as hdu:
    bias = CCDData(hdu[0].data, unit=u.adu)

# Loop through the image file collection
for image, fname in ifc.hdus(return_fname=True):
    data = CCDData(image.data, unit=u.adu)

    # Bias Subtract
    bias_subtracted = subtract_bias(data, bias)

    # Trim
    trimmed = bias_subtracted[:, np.r_[xdata_slice_list]]
    trimmed = trimmed[np.r_[ydata_slice_list], :]
    
    #Cutout
    y_center = trimmed.data.shape[0] / 2
    x_center = trimmed.data.shape[1] / 2
    position = (x_center, y_center)
#   size = (2640, 2640)
    size = (10560, 10560)
    cutout = Cutout2D(trimmed.data, position, size)
    cutout = CCDData(cutout.data, unit=u.adu)
    
    # reduce filesize float64 -> float32
    cutout.data = cutout.data.astype('float32')

    # Save new image
    filename = './trimmed/' + 'trimmed_' + str(fname)
    cutout.write(filename)
