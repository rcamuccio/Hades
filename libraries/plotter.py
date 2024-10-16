from config import Configuration

import matplotlib.pyplot as plt

class Plotter:

	@staticmethod
	def plot_airmass(x_list, y_list):

		plt.clf()
		font = {"fontname":"Monospace", "size":12}

		plt.figure(figsize=(10,10))
		plt.plot(x_list, y_list, color="black")

		plt.title("Air Mass", **font)
		plt.xlabel("Time [JD]", **font)
		plt.ylabel("Air Mass", **font)
		plt.xticks(**font)
		plt.yticks(**font)

		plt.savefig("plot-airmass.png", dpi=400)
		plt.close()

		return

	@staticmethod
	def plot_colormag(name, x_list, y_list, yfit):

		plt.clf()
		font = {"fontname":"Monospace", "size":12}

		plt.figure(figsize=(10,10))
		plt.scatter(x_list, y_list, s=1, color="gray")
		plt.plot(x_list, yfit, color="blue")

		plt.title("Color-Magnitude Diagram", **font)
		plt.xlabel("Color index [mag]", **font)
		plt.ylabel("Delta magnitude [mag]", **font)
		plt.xticks(**font)
		plt.yticks(**font)

		plt.savefig(name + "-colormag.png", dpi=400)
		plt.close()

		return

	@staticmethod
	def plot_extinction(name, x_list, y_list, yfit):

		plt.clf()
		font = {"fontname":"Monospace", "size":12}

		plt.figure(figsize=(10,10))
		plt.scatter(x_list, y_list, s=1, color="gray")
		plt.plot(x_list, yfit, color="blue")

		plt.title("Extinction", **font)
		plt.xlabel("Air mass", **font)
		plt.ylabel("Instrumental magnitude [mag]", **font)
		plt.xticks(**font)
		plt.yticks(**font)

		plt.savefig(name + "-extinction.png", dpi=400)
		plt.close()

		return

	@staticmethod
	def plot_field(name, data, apertures, boxes, save=True):

		#interval = MinMaxInterval()
		#vmin, vmax = interval.get_limits(data)
		#norm = ImageNormalize(vmin=vmin/-134, vmax=vmax/263, stretch=LinearStretch())

		plt.figure(figsize=(10,10))
		plt.imshow(data, cmap="inferno", origin="lower", norm=norm, interpolation="nearest")

		apertures.plot(color="lime", lw=0.5, alpha=0.5)

		mask_box = boxes["mask_box"]
		eml_box = boxes["eml_box"]
		emt_box = boxes["emt_box"]
		emr_box = boxes["emr_box"]
		emb_box = boxes["emb_box"]

		mask_box.plot(color="red", ls="dashed")
		eml_box.plot(color="cyan", ls="dashed")
		emt_box.plot(color="lime", ls="dashed")
		emr_box.plot(color="yellow", ls="dashed")
		emb_box.plot(color="orange", ls="dashed")

		if save:
			plt.savefig(name + "-field.png", dpi=400)
		else:
			plt.show()
		plt.close()

		return

	@staticmethod
	def plot_growth_radius(x_list, y_list):

		plt.clf()
		font = {"fontname":"Monospace", "size":12}

		plt.figure(figsize=(10,10))
		plt.plot(x_list, y_list, color="black")

		plt.title("Growth Radius", **font)
		plt.xlabel("Time [JD]", **font)
		plt.ylabel("Radius [px]", **font)
		plt.xticks(**font)
		plt.yticks(**font)

		plt.savefig("plot-growth-radius.png", dpi=400)
		plt.close()

		return

	@staticmethod
	def plot_histogram(name, x_list):

		plt.clf()
		font = {"fontname":"Monospace", "size":12}

		plt.figure(figsize=(10,10))
		plt.hist(x_list, bins=100)

		plt.savefig(name + "-histogram.png", dpi=400)
		plt.close()

		return

	@staticmethod
	def plot_lightcurve(name, x_list, y_list, err_list):

		plt.clf()
		font = {"fontname":"Monospace", "size":12}

		plt.figure(figsize=(10,10))
		plt.errorbar(x_list, y_list, yerr=err_list, fmt="o", linewidth=0.5, markersize=0.5, capsize=2, capthick=0.5)

		#plt.ylim(13, 19)
		plt.title("Time Series of " + name, **font)
		plt.xlabel("Time [JD]", **font)
		plt.ylabel("Magnitude [mag]", **font)
		plt.xticks(**font)
		plt.yticks(**font)

		plt.savefig(name + "-timeseries.png", dpi=400)
		plt.close()

		return

	@staticmethod
	def plot_magerr(name, x_list, y_list):

		plt.clf()
		font = {"fontname":"Monospace", "size":12}

		plt.figure(figsize=(10,10))
		plt.scatter(x_list, y_list, s=1, color="gray")

		plt.title("Magnitude-Error Diagram", **font)
		plt.xlabel("Instrumental Magnitude [mag]", **font)
		plt.ylabel("Instrumental Magnitude Error [mag]", **font)
		plt.xticks(**font)
		plt.yticks(**font)

		plt.savefig(name + "-magerror.png", dpi=400)
		plt.close()

		return

	@staticmethod
	def plot_seeing_pix(x_list, y_list):

		plt.clf()
		font = {"fontname":"Monospace", "size":12}

		plt.figure(figsize=(10,10))
		plt.plot(x_list, y_list, color="black")

		plt.title("Seeing (Pixels)", **font)
		plt.xlabel("Time [JD]", **font)
		plt.ylabel("Seeing [px]", **font)
		plt.xticks(**font)
		plt.yticks(**font)

		plt.savefig("plot-seeing-pix.png", dpi=400)
		plt.close()

		return

	@staticmethod
	def plot_seeing_sky(x_list, y_list):

		plt.clf()
		font = {"fontname":"Monospace", "size":12}

		plt.figure(figsize=(10,10))
		plt.plot(x_list, y_list, color="black")

		plt.title("Seeing (Sky)", **font)
		plt.xlabel("Time [JD]", **font)
		plt.ylabel("Seeing [arcsec]", **font)
		plt.xticks(**font)
		plt.yticks(**font)

		plt.savefig("plot-seeing-sky.png", dpi=400)
		plt.close()

		return