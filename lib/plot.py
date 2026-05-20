from config import Configuration
from astropy.stats import sigma_clipped_stats
from astropy.visualization import ZScaleInterval
import matplotlib.pyplot as plt

class Plot:

	@staticmethod
	def airmass(xlist, ylist):

		plt.clf()
		font = {'fontname':'Monospace', 'size':12}

		plt.figure(figsize=(10, 10))
		plt.plot(xlist, ylist, color='black')

		plt.title('Air Mass', **font)
		plt.xlabel('Time [JD]', **font)
		plt.ylabel('Air Mass', **font)
		plt.xticks(**font)
		plt.yticks(**font)

		plt.savefig('plot-airmass.png', dpi=400)
		plt.close()

	@staticmethod
	def colormag(xlist, ylist, yfit):

		plt.clf()
		font = {'fontname':'Monospace', 'size':12}

		plt.figure(figsize=(10, 10))
		plt.scatter(xlist, ylist, s=1, color='gray')
		plt.plot(xlist, yfit, color='blue')

		plt.title('Color-Magnitude Diagram', **font)
		plt.xlabel('Color index [mag]', **font)
		plt.ylabel('Delta magnitude [mag]', **font)
		plt.xticks(**font)
		plt.yticks(**font)

		plt.savefig('plot-colormag.png', dpi=400)
		plt.close()

	@staticmethod
	def extinction(xlist, ylist, yfit):

		plt.clf()
		font = {'fontname':'Monospace', 'size':12}

		plt.figure(figsize=(10, 10))
		plt.scatter(xlist, ylist, s=1, color='gray')
		plt.plot(xlist, yfit, color='blue')

		plt.title('Extinction', **font)
		plt.xlabel('Air mass', **font)
		plt.ylabel('Instrumental magnitude [mag]', **font)
		plt.xticks(**font)
		plt.yticks(**font)

		plt.savefig('plot-extinction.png', dpi=400)
		plt.close()

	@staticmethod
	def field(data, path, apertures=None, boxes=None, save=Configuration.SAVE_FIGURE):

		plt.clf()
		font = {'fontname':Configuration.FONT_NAME, 'size':Configuration.FONT_SIZE}
		plt.figure(figsize=Configuration.FIGURE_SIZE)
		interval = ZScaleInterval()
		vmin, vmax = interval.get_limits(data)
		plt.imshow(data, cmap=Configuration.CMAP, origin='lower', vmin=vmin, vmax=vmax)
		plt.colorbar()
		plt.xlabel('Pixel Number [x]')
		plt.ylabel('Pixel Number [y]')

		if apertures != None:
			apertures.plot(color='lime', lw=0.5, alpha=0.5)

		if boxes != None:
			eml_box = boxes['eml_box']
			emt_box = boxes['emt_box']
			emr_box = boxes['emr_box']
			emb_box = boxes['emb_box']

			eml_box.plot(color='cyan', ls='dashed')
			emt_box.plot(color='cyan', ls='dashed')
			emr_box.plot(color='cyan', ls='dashed')
			emb_box.plot(color='cyan', ls='dashed')

		if save:
			plt.savefig(path, dpi=Configuration.DPI)

		plt.close()

	@staticmethod
	def growth_radius(xlist, ylist):

		plt.clf()
		font = {'fontname':'Monospace', 'size':12}

		plt.figure(figsize=(10, 10))
		plt.plot(xlist, ylist, color='black')

		plt.title('Growth Radius', **font)
		plt.xlabel('Time [JD]', **font)
		plt.ylabel('Radius [px]', **font)
		plt.xticks(**font)
		plt.yticks(**font)

		plt.savefig('plot-growth-radius.png', dpi=400)
		plt.close()

	@staticmethod
	def histogram(data, path, save=Configuration.SAVE_FIGURE):

		plt.clf()
		font = {'fontname':Configuration.FONT_NAME, 'size':Configuration.FONT_SIZE}
		plt.figure(figsize=Configuration.FIGURE_SIZE)

		# calculate the histogram limits
		data_mn, data_md, data_sd = sigma_clipped_stats(data, sigma=Configuration.SIG_BKG)
		xmin = data_md - Configuration.HISTOGRAM_LIMIT * data_sd
		xmax = data_md + Configuration.HISTOGRAM_LIMIT * data_sd

		plt.hist(data.flatten(), bins=Configuration.HISTOGRAM_BINS, range=(xmin, xmax), histtype=Configuration.HISTOGRAM_TYPE)

		plt.yscale(Configuration.HISTOGRAM_SCALE)

		plt.xlabel('Pixel Value [ADU]')
		plt.ylabel('Pixel Count')

		plt.savefig(path, dpi=400)
		plt.close()

	@staticmethod
	def lightcurve(xlist, ylist, elist):

		plt.clf()
		font = {'fontname':'Monospace', 'size':12}

		plt.figure(figsize=(10, 10))
		plt.errorbar(xlist, ylist, yerr=elist, fmt='o', linewidth=0.5, markersize=0.5, capsize=2, capthick=0.5)

		#plt.ylim(13, 19)
		plt.title('Time Series of ' + name, **font)
		plt.xlabel('Time [JD]', **font)
		plt.ylabel('Magnitude [mag]', **font)
		plt.xticks(**font)
		plt.yticks(**font)

		plt.savefig('plot-' + name + '-timeseries.png', dpi=400)
		plt.close()

	@staticmethod
	def seeing(xlist, ylist, type):

		plt.clf()
		font = {'fontname':'Monospace', 'size':12}

		plt.figure(figsize=(10, 10))
		plt.plot(xlist, ylist, color='black')

		if type == 'physical':
			plt.title('Seeing (physical)', **font)
			plt.xlabel('Time [JD]', **font)
			plt.ylabel('Seeing [px]', **font)
		elif type == 'sky':
			plt.title('Seeing (sky)', **font)
			plt.xlabel('Time [JD]', **font)
			plt.ylabel('Seeing [arcsec]', **font)

		plt.xticks(**font)
		plt.yticks(**font)
		
		plt.savefig('plot-seeing.png', dpi=400)
		plt.close()