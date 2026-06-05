from config import Configuration

import numpy as np

class Observatory:

	@staticmethod
	def axis_x():

		if Configuration.OBSERVATORY == 'toros':
			axis_x = 10560
		else:
			axis_x = None

		return axis_x

	@staticmethod
	def axis_x_raw():

		if Configuration.OBSERVATORY == 'toros':
			axis_x_raw = 12000
		else:
			axis_x_raw = None

		return axis_x_raw

	@staticmethod
	def axis_y():

		if Configuration.OBSERVATORY == 'toros':
			axis_y = 10560
		else:
			axis_y = None

		return axis_y

	@staticmethod
	def axis_y_raw():

		if Configuration.OBSERVATORY == 'toros':
			axis_y_raw = 10600
		else:
			axis_y_raw = None

		return axis_y_raw

	@staticmethod
	def bandpass():

		if Configuration.OBSERVATORY == 'toros':
			bandpass = [147, 141, 147, 147]
		else:
			bandpass = None

		return bandpass

	@staticmethod
	def binning():

		if Configuration.OBSERVATORY == 'toros':
			binning = 1

		return binning

	@staticmethod
	def center_wavelength():

		if Configuration.OBSERVATORY == 'toros':
			center_wavelength = [473.5, 638.5, 775.5, 922.5]
		else:
			center_wavelength = None

		return center_wavelength

	@staticmethod
	def declination_limit():

		if Configuration.OBSERVATORY == 'toros':
			declination_limit = 26.66
		else:
			declination_limit = None

		return declination_limit

	@staticmethod
	def elevation():

		if Configuration.OBSERVATORY == 'toros':
			elevation = 2420
		else:
			elevation = None

		return elevation

	@staticmethod
	def exposure_number():

		if Configuration.OBSERVATORY == 'toros':
			exposure_number = 1
		else:
			exposure_number = None

		return exposure_number

	@staticmethod
	def exposure_time(type='light', scale='s'):

		if Configuration.OBSERVATORY == 'toros':
			if type == 'dark-f':
				exposure_time = 5.
			elif type == 'dark-l':
				exposure_time = 300.
			elif type == 'flat':
				exposure_time = 5.
			elif type == 'light':
				exposure_time = 300.

			if scale == 'd':
				exposure_time = exposure_time / 60. / 60. / 24.

		return exposure_time

	@staticmethod
	def fov():

		if Configuration.OBSERVATORY == 'toros':
			pixel_scale = Observatory.pixel_scale()
			axis_x = Observatory.axis_x()
			fov = (pixel_scale * axis_x) / 3600.
		else:
			fov = None

		return fov

	@staticmethod
	def gain():

		if Configuration.OBSERVATORY == 'toros':
			gain = 0.380
		else:
			gain = None

		return gain

	@staticmethod
	def latitude():

		if Configuration.OBSERVATORY == 'toros':
			latitude = -31.8023
		else:
			latitude = None

		return latitude

	@staticmethod
	def longitude():

		if Configuration.OBSERVATORY == 'toros':
			longitude = -69.3265
		else:
			longitude = None

		return longitude

	@staticmethod
	def mirror_diameter():

		if Configuration.OBSERVATORY == 'toros':
			mirror_diameter = 0.610
		else:
			mirror_diameter = None

		return mirror_diameter

	@staticmethod
	def mirror_radius():

		if Configuration.OBSERVATORY == 'toros':
			mirror_diameter = Observatory.mirror_diameter()
			mirror_radius = mirror_diameter / 2
		else:
			mirror_radius = None

		return mirror_radius

	@staticmethod
	def overhead():

		if Configuration.OBSERVATORY == 'toros':
			overhead = 30.
		else:
			overhead = None

		return overhead

	@staticmethod
	def overscan_x():

		if Configuration.OBSERVATORY == 'toros':
			overscan_x = 180
		else:
			overscan_x = None

		return overscan_x

	@staticmethod
	def overscan_y():

		if Configuration.OBSERVATORY == 'toros':
			overscan_y = 20
		else:
			overscan_y = None

		return overscan_y

	@staticmethod
	def peakmax():

		if Configuration.OBSERVATORY == 'toros':
			peakmax = 45000.
		else:
			peakmax = None

		return peakmax

	@staticmethod
	def pixel_scale():

		if Configuration.OBSERVATORY == 'toros':
			pixel_scale = 0.4959
		else:
			pixel_scale = None

		return pixel_scale

	@staticmethod
	def quantum_efficiency_atm():

		if Configuration.OBSERVATORY == 'toros':
			quantum_efficiency_atm = [0.8, 0.9, 0.9, 0.9]
		else:
			quantum_efficiency_atm = None

		return quantum_efficiency_atm

	@staticmethod
	def quantum_efficiency_ccd():

		if Configuration.OBSERVATORY == 'toros':
			quantum_efficiency_ccd = 0.85
		else:
			quantum_efficiency_ccd = None

		return quantum_efficiency_ccd

	@staticmethod
	def quantum_efficiency_fil():

		if Configuration.OBSERVATORY == 'toros':
			quantum_efficiency_fil = 0.9
		else:
			quantum_efficiency_fil = None

		return quantum_efficiency_fil

	@staticmethod
	def quantum_efficiency_pri():

		if Configuration.OBSERVATORY == 'toros':
			quantum_efficiency_pri = 0.96
		else:
			quantum_efficiency_pri = None

		return quantum_efficiency_pri

	@staticmethod
	def quantum_efficiency_sec():

		if Configuration.OBSERVATORY == 'toros':
			quantum_efficiency_sec = 0.96
		else:
			quantum_efficiency_sec = None

		return quantum_efficiency_sec

	@staticmethod
	def read_noise():

		if Configuration.OBSERVATORY == 'toros':
			read_noise = 5.
		else:
			read_noise = None

		return read_noise

	@staticmethod
	def read_time():

		if Configuration.OBSERVATORY == 'toros':
			read_time = 90.
		else:
			read_time = None

		return read_time

	@staticmethod
	def seeing():

		if Configuration.OBSERVATORY == 'toros':
			seeing = 0.93
		else:
			seeing = None

		return seeing

	@staticmethod
	def sky():

		if Configuration.OBSERVATORY == 'toros':
			sky = [22.1, 21.1, 20.1, 18.7]
		else:
			sky = None

		return sky

	@staticmethod
	def tap_x():

		if Configuration.OBSERVATORY == 'toros':
			tap_x = 1320
		else:
			tap_x = None

		return tap_x

	@staticmethod
	def tap_y():

		if Configuration.OBSERVATORY == 'toros':
			tap_y = 5280
		else:
			tap_y = None

		return tap_y

	@staticmethod
	def throughput():

		if Configuration.OBSERVATORY == 'toros':
			qe_atm = Observatory.quantum_efficiency_atm()
			qe_atm = float(np.mean(qe_atm))
			qe_ccd = Observatory.quantum_efficiency_ccd()
			qe_fil = Observatory.quantum_efficiency_fil()
			qe_pri = Observatory.quantum_efficiency_pri()
			qe_sec = Observatory.quantum_efficiency_sec()
			vignetting = Observatory.vignetting()

			throughput = qe_atm * qe_ccd * qe_fil * qe_pri * qe_sec * vignetting

		else:
			throughput = None

		return throughput

	@staticmethod
	def total_exposure():

		if Configuration.OBSERVATORY == 'toros':
			exposure_number = Observatory.exposure_number()
			exposure_time = Observatory.exposure_time()
			overhead = Observatory.overhead()
			read_time = Observatory.read_time()
			total_exposure = (exposure_time + read_time) * exposure_number + overhead
		else:
			total_exposure = None

		return total_exposure

	@staticmethod
	def utc():

		if Configuration.OBSERVATORY == 'toros':
			utc = -3
		else:
			utc = None

		return utc

	@staticmethod
	def vignetting():

		if Configuration.OBSERVATORY == 'toros':
			vignetting = 0.756
		else:
			vignetting = None

		return vignetting