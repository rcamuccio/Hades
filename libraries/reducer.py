from config import Configuration

from astropy import units as u
from astropy.io import fits
from astropy.stats import SigmaClip, sigma_clipped_stats
from astropy.wcs import WCS
from astropy.wcs.utils import pixel_to_skycoord
from photutils.aperture import RectangularAperture
from photutils.background import Background2D, MedianBackground
from photutils.detection import DAOStarFinder
from photutils.segmentation import detect_threshold, detect_sources
from photutils.utils import circular_footprint
from reproject import reproject_interp
import ccdproc
import glob
import numpy as np
import os
import shutil
import subprocess

class Reducer:

	@staticmethod
	def align_frames(obj_dir):

		wcs_dir = os.path.join(obj_dir, 'wcs')
		align_dir = os.path.join(obj_dir, 'align')

		os.chdir(wcs_dir)
		obj_list = []
		for item in glob.glob('*.fit'):
			obj_list.append(item)
		obj_list = sorted(obj_list)

		reference_frame = fits.open(obj_list[0])
		reference_data = reference_frame[0].data
		reference_header = reference_frame[0].header

		align_path = os.path.join(align_dir, 'a-' + obj_list[0])

		if not os.path.isfile(align_path):
			print('Saving reference frame', obj_list[0])
			reference_hdu = fits.PrimaryHDU(reference_data, header=reference_header)
			reference_hdu.writeto(align_path, overwrite=True)

		for i in range(1, len(obj_list)):

			align_path = os.path.join(align_dir, 'a-' + obj_list[i])

			if os.path.isfile(align_path):
				print('Skipping alignment on frame', obj_list[i])

			else:
				print('Aligning frame', obj_list[i], 'with reference frame')
				target_frame = fits.open(obj_list[i])
				target_data = target_frame[0].data
				target_header = target_frame[0].header

				target_hdu = fits.PrimaryHDU(target_data, header=target_header)
				array, footprint = reproject_interp(target_hdu, reference_header)

				aligned_hdu = fits.PrimaryHDU(array, header=reference_header)
				aligned_hdu.writeto(align_path)

		return obj_list

	@staticmethod
	def make_dark(dark_dir):

		dark_path = os.path.join(dark_dir, 'master-dark.fit')

		if os.path.isfile(dark_path):
			print('Reading extant master dark')
			master_dark = fits.open(dark_path)

		else:
			print('Creating master dark')
			os.chdir(dark_dir)
			
			dark_list = []

			for item in glob.glob('*.fit'):
				dark_list.append(item)

			master_dark = ccdproc.combine(dark_list, method=Configuration.COMBINE_METHOD, unit='adu', mem_limit=Configuration.MEM_LIMIT, dtype=Configuration.DTYPE)

			ccdproc.fits_ccddata_writer(master_dark, dark_path, overwrite=True)

		return master_dark

	@staticmethod
	def make_flat(flat_dir, dark_dir):

		flat_path = os.path.join(flat_dir, 'flatfield.fit')
		dark_path = os.path.join(dark_dir, 'master-dark.fit')

		if os.path.isfile(flat_path):
			print('Reading extant flatfield')
			flatfield = fits.open(flat_path)

		else:
			print('Creating flatfield')
			master_dark = ccdproc.fits_ccddata_reader(dark_path)
			dark_exposure = master_dark.header['exposure']*u.second

			os.chdir(flat_dir)

			flat_list = []
			for item in glob.glob('*.fit'):

				flat = ccdproc.fits_ccddata_reader(item, unit='adu')
				flat_exposure = flat.header['exposure']*u.second

				reduced_flat = ccdproc.subtract_dark(flat, master_dark, data_exposure=flat_exposure, dark_exposure=dark_exposure)
				flat_list.append(reduced_flat)

			combined_flat = ccdproc.combine(flat_list, method=Configuration.COMBINE_METHOD, unit='adu', mem_limit=Configuration.MEM_LIMIT, dtype=Configuration.DTYPE)

			combined_flat_data = np.asarray(combined_flat)
			flatfield_data = combined_flat_data / np.mean(combined_flat_data)

			flatfield = ccdproc.CCDData(flatfield_data, unit='adu')
			ccdproc.fits_ccddata_writer(flatfield, flat_path)

		return flatfield

	@staticmethod
	def make_mask(object_frame):

		frame = fits.open(object_frame)
		frame_data = frame[0].data
		frame_header = frame[0].header

		mask = np.zeros(frame_data.shape, dtype=bool)

		binning = frame_header['XBINNING']
		edge = 100 // binning

		dx = frame_data.shape[0]
		dy = frame_data.shape[1]

		# --- Bad pixel mask
		mask_y1 = 860 // binning
		mask_y2 = 4089 // binning
		mask_delta_y = mask_y2 - mask_y1
		center_y = mask_delta_y // 2

		mask_x1 = 1200 // binning
		mask_x2 = 1210 // binning
		mask_delta_x = mask_x2 - mask_x1
		center_x = mask_delta_x // 2

		mask[mask_y1:mask_y2, mask_x1:mask_x2] = True

		# --- Left edge mask
		eml_center_x = edge // 2
		eml_center_y = dy // 2
		mask[0:dy, 0:edge] = True

		# --- Top edge mask
		emt_center_x = dx // 2
		emt_center_y = edge // 2
		mask[dy-edge:dy, 0:dx] = True

		# --- Right edge mask
		emr_center_x = edge // 2
		emr_center_y = dy // 2
		mask[0:dy, dx-edge:dx] = True

		# --- Bottom edge mask
		emb_center_x = dx // 2
		emb_center_y = edge // 2
		mask[0:edge, 0:dx] = True

		boxes = {}

		mask_box = RectangularAperture((center_x + mask_x1, center_y + mask_y1), mask_delta_x, mask_delta_y, theta=0.)
		eml_box = RectangularAperture((eml_center_x, eml_center_y), edge, dy, theta=0.)
		emt_box = RectangularAperture((emt_center_x, emt_center_y + dy - edge), dx, edge, theta=0.)
		emr_box = RectangularAperture((emr_center_x + dx - edge, emr_center_y), edge, dy, theta=0.)
		emb_box = RectangularAperture((emb_center_x, emb_center_y), dx, edge, theta=0.)

		boxes['mask_box'] = mask_box
		boxes['eml_box'] = eml_box
		boxes['emt_box'] = emt_box
		boxes['emr_box'] = emr_box
		boxes['emb_box'] = emb_box

		return mask, boxes

	@staticmethod
	def make_stack(obj_dir):

		align_dir = os.path.join(obj_dir, 'align')
		stack_path = os.path.join(align_dir, 'stack.fit')

		if os.path.isfile(stack_path):
			print('Reading extant stack')
			stack = ccdproc.fits_ccddata_reader(stack_path)

		else:
			print('Creating master stack')
			os.chdir(align_dir)
			stack_list = []
			for item in glob.glob('*.fit'):
				stack_list.append(item)

			stack = ccdproc.combine(stack_list, method=Configuration.COMBINE_METHOD, unit='adu', mem_limit=Configuration.MEM_LIMIT, dtype=Configuration.DTYPE)
			ccdproc.fits_ccddata_writer(stack, stack_path)

		return stack

	@staticmethod
	def reduce_objects(obj_dir, flat_dir, dark_dir, bkg_method='flat'):

		flat_path = os.path.join(flat_dir, 'flatfield.fit')
		dark_path = os.path.join(dark_dir, 'master-dark.fit')

		raw_dir = os.path.join(obj_dir, 'raw')
		cal_dir = os.path.join(obj_dir, 'cal')

		os.chdir(raw_dir)
		obj_list = []
		for item in glob.glob('*.fit'):
			obj_list.append(item)
		obj_list = sorted(obj_list)

		for obj in obj_list:

			obj_name = 'red-' + obj
			cal_path = os.path.join(cal_dir, obj_name)

			if os.path.isfile(cal_path):
				print('Skipping reduction on frame', obj)

			else:
				print('Reducing frame', obj)
				obj_frame = fits.open(obj)
				obj_frame_data = obj_frame[0].data
				obj_frame_header = obj_frame[0].header

				master_dark = fits.open(dark_path)
				master_dark_data = master_dark[0].data
				master_dark_header = master_dark[0].header
				reduced_obj_frame_data = obj_frame_data - master_dark_data

				flatfield = fits.open(flat_path)
				flatfield_data = flatfield[0].data
				flatfield_header = flatfield[0].header
				reduced_obj_frame_data /= flatfield_data

				sigma_clip = SigmaClip(sigma=Configuration.SIGMA)
				threshold = detect_threshold(reduced_obj_frame_data, nsigma=Configuration.SIGMA_BKG, sigma_clip=sigma_clip)
				segment_img = detect_sources(reduced_obj_frame_data, threshold, npixels=Configuration.NPIXELS)
				footprint = circular_footprint(radius=10)
				mask = segment_img.make_source_mask(footprint=footprint)

				bkg_estimator = MedianBackground()
				bkg = Background2D(reduced_obj_frame_data, box_size=Configuration.BOX_SIZE, filter_size=Configuration.FILTER_SIZE, sigma_clip=sigma_clip, bkg_estimator=bkg_estimator)
				
				print(bkg.background_median, bkg.background_rms_median)

				if bkg_method == '2d':
					reduced_obj_frame_data -= bkg.background

				elif bkg_method == 'flat':
					reduced_obj_frame_data -= bkg.background_median

				else:
					reduced_obj_frame_data -= bkg.background

				obj_hdu = fits.PrimaryHDU(reduced_obj_frame_data, header=obj_frame_header)
				obj_hdu.writeto(cal_path)

		return obj_list

	@staticmethod
	def solve_plates(obj_dir):

		cal_dir = os.path.join(obj_dir, 'cal')
		wcs_dir = os.path.join(obj_dir, 'wcs')

		os.chdir(cal_dir)
		obj_list = []
		for item in glob.glob('*.fit'):
			obj_list.append(item)
		obj_list = sorted(obj_list)

		for obj in obj_list:

			cal_path = os.path.join(cal_dir, obj)
			wcs_path = os.path.join(wcs_dir, 'wcs-' + obj)

			if os.path.isfile(wcs_path):
				print('Skipping plate solve on frame', obj)

			else:
				print('Plate solving frame', obj)
				file_name = obj[:-4]
				file_axy = file_name + '.axy'
				file_corr = file_name + '.corr'
				file_match = file_name + '.match'
				file_new = file_name + '.new'
				file_rdls = file_name + '.rdls'
				file_solved = file_name + '.solved'
				file_wcs = file_name + '.wcs'
				file_xyls = file_name + '-indx.xyls'

				subprocess.run(['solve-field', '--no-plots', obj, '--ra', str(Configuration.FIELD_RA), '--dec', str(Configuration.FIELD_DEC), '--radius', str(Configuration.RAD_SOLVE)])
				subprocess.run(['rm', file_axy])
				subprocess.run(['rm', file_corr])
				subprocess.run(['rm', file_match])
				subprocess.run(['rm', file_rdls])
				subprocess.run(['rm', file_solved])
				subprocess.run(['rm', file_wcs])
				subprocess.run(['rm', file_xyls])
				subprocess.run(['mv', file_new, 'wcs-' + str(file_name) + '.fit'])

		wcs_list = []
		for item in glob.glob('wcs-*.fit'):

			old_item_path = os.path.join(cal_dir, item)
			new_item_path = os.path.join(wcs_dir, item)

			shutil.move(old_item_path, new_item_path)

		return obj_list
