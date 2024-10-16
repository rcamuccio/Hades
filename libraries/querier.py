from config import Configuration

from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import ascii
from astropy.table import Table
from astroquery.gaia import Gaia
import numpy as np
import requests

class Querier:

	@staticmethod
	def ps1_checklegal(table, release):

		releaselist = ("dr1", "dr2")

		if release not in ("dr1", "dr2"):
			raise ValueError("Bad value for release (must be one of {})".format(", ".join(releaselist)))

		if release == "dr1":
			tablelist = ("mean", "stack")
		else:
			tablelist = ("mean", "stack", "detection")

		if table not in tablelist:
			raise ValueError("Bad value for table (for {} must be one of {})".format(release, ", ".join(tablelist)))

	@staticmethod
	def ps1_cone(ra, dec, radius, table="mean", release="dr2", format="csv", columns=None, baseurl="https://catalogs.mast.stsci.edu/api/v0.1/panstarrs", verbose=False, **kw):

		data = kw.copy()
		data["ra"] = ra
		data["dec"] = dec
		data["radius"] = radius

		return Querier.ps1_search(table=table, release=release, format=format, columns=columns, baseurl=baseurl, verbose=verbose, **data)

	@staticmethod
	def ps1_metadata(table="mean", release="dr2", baseurl="https://catalogs.mast.stsci.edu/api/v0.1/panstarrs"):

		Querier.ps1_checklegal(table, release)

		url = f"{baseurl}/{release}/{table}/metadata"

		r = requests.get(url)
		r.raise_for_status()
		v = r.json()

		table = Table(rows=[(x["name"], x["type"], x["description"]) for x in v], names=("name", "type", "description"))

		return table

	@staticmethod
	def ps1_search(table="mean", release="dr2", format="csv", columns=None, baseurl="https://catalogs.mast.stsci.edu/api/v0.1/panstarrs", verbose=False, **kw):

		data = kw.copy()

		if not data:
			raise ValueError("You must specify some parameters for search")

		Querier.ps1_checklegal(table, release)

		if format not in ("csv", "votable", "json"):
			raise ValueError("Bad value for format")

		url = f"{baseurl}/{release}/{table}.{format}"

		if columns:
			dcols = {}

			for col in Querier.ps1_metadata(table, release)["name"]:
				dcols[col.lower()] = 1

			badcols = []
			for col in columns:
				if col.lower().strip() not in dcols:
					badcols.append(col)
			if badcols:
				raise ValueError("Some columns not found in table: {}".format(", ".join(badcols)))

			data["columns"] = "[{}]".format(",".join(columns))

		r = requests.get(url, params=data)

		if verbose:
			print(r.url)

		r.raise_for_status()

		if format == "json":
			return r.json()
		else:
			return r.text

	@staticmethod
	def query_gaia_cone(ra, dec, radius):

		print("Submitting Gaia cone search")
		Gaia.MAIN_GAIA_TABLE = "gaiadr3.gaia_source"
		Gaia.ROW_LIMIT = -1

		ra = str(ra)
		dec = str(dec)

		search_coord = SkyCoord(ra, dec, unit=(u.deg, u.deg), frame="fk5")
		search_radius = u.Quantity(radius, u.deg)

		print(search_coord)
		print(search_radius)

		search = Gaia.cone_search_async(search_coord, radius=search_radius)

		table = search.get_results()

		return table

	@staticmethod
	def query_gaia_square(ra, dec, side):

		print("Submitting Gaia square search")
		Gaia.MAIN_GAIA_TABLE = "gaiadr2.gaia_source"
		Gaia.ROW_LIMIT = -1

		ra = str(ra)
		dec = str(dec)
		side = float(side)

		search_coord = SkyCoord(ra, dec, unit=(u.deg, u.deg), frame="fk5")
		search_width = u.Quantity(side, u.deg)
		search_height = u.Quantity(side, u.deg)

		table = Gaia.query_object_async(search_coord, width=search_width, height=search_height)

		return table

	@staticmethod
	def query_ps1_cone(ra, dec, radius):

		print("Submitting PS1 cone search")
		release = "dr2"
		constraints = {"nDetections.gt":1}

		ra = float(ra)
		dec = float(dec)
		radius = float(radius)

		columns = """objID, raMean, decMean, nDetections, ng, nr, ni, nz, ny, gMeanPSFMag, rMeanPSFMag, iMeanPSFMag, zMeanPSFMag, yMeanPSFMag""".split(",")
		columns = [x.strip() for x in columns]
		columns = [x for x in columns if x and not x.startswith("#")]

		results = Querier.ps1_cone(ra, dec, radius, release=release, columns=columns, **constraints)

		table = ascii.read(results)

		for filter in "grizy":
			col = filter + "MeanPSFMag"
			try:
				table[col].format = ".4f"
				table[col][table[col] == -999.0] = np.nan
			except KeyError:
				print("{} not found".format(col))	

		return table