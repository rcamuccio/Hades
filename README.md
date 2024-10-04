# Hades

**HADES** is an astronomical observer's toolbox. It contains a set of tools designed to accomplish a variety of observational tasks, including image processing, aperture photometry, and responding to astronomical alerts.

---

## Usage

The code is run through a set of scripts from the main HADES directory. The scripts currently available include:

- GCN listener (`gcn_listener.py`)
- Limiting magnitude (`limiting_magnitude.py`)
- Quick reduction (`quick_reduction.py`)
- Relative photometry (`relative_photometry.py`)

### Configuration

The settings are controlled in the `config.py` file.

### Example 1: Running the GCN listener

```
$ python -m scripts.gcn_listener
```

### Example 2: Running quick reduction on observational data

```
$ python -m scripts.quick_reduction
```

The quick reduction script is run on a single night (yyyy-mm-dd) of data, assuming the following directory format:

```
[yyyy-mm-dd]
	[dark]
		[dark-flat]
			image1.fit
			image2.fit
			...
		[dark-obj]
			image1.fit
			image2.fit
			...
	[flat]
		image1.fit
		image2.fit
		...
	[obj]
		image1.fit
		image2.fit
		...
```

---

## Installation

### Development Platform

HADES was developed on a machine (`Epimetheus`) running Ubuntu 22.04.4 LTS and using an Anaconda environment with Python 3.12.3. The list of library requirements are listed in the the `requirements.txt` file in the main directory.

### GCN Subscription

One requires an account on the NASA General Coordinates Network (GCN) platform (https://gcn.nasa.gov). The parameters `client_id` and `client_secret` are uniquely generated per user, and must be plugged into the configuration file in order to use the GCN listener.

### Catalogs

The GCN listener uses the GLADE catalog (both the latest 'GLADE+' and previous 'GLADE 2.4' catalogs, https://glade.elte.hu/) for selecting galaxy targets to follow up. Download the catalogs as ASCII text files (e.g. GLADE+.txt and/or GLADE_2.4.txt) and save them in a separate directory.

### Astrometry.net

The reduction toolbox includes plate solving FITS files. The plate solving routine is handled by a local instance of the Astrometry.net (http://astrometry.net/use.html) code. One needs to properly install Astrometry.net and the index files it uses for solving fields. There is a shell script available (https://github.com/rcamuccio/package-installers/blob/master/astrometry-installer.sh) which can be used to download all of the index files required. Please note that the shell script is still in development, and will likely need to be modified depending on the user's configuration.

---

## References

- G. Dálya et al. (2018), GLADE: A Galaxy Catalogue for Multimessenger Searches in the Advanced Gravitational-Wave Detector Era, MNRAS, Volume 479, Issue 2, https://doi.org/10.1093/mnras/sty1703

- G. Dálya et al. (2022), GLADE+: An Extended Galaxy Catalogue for Multimessenger Searches with Advanced Gravitational-Wave Detectors, MNRAS, Volume 514, Issue 1, https://doi.org/10.1093/mnras/stac1443

- LVK Collaboration (2022), LIGO/Virgo/KAGRA Public Alerts User Guide, https://emfollow.docs.ligo.org/userguide/

---

7 Jun 2019<br>
Last update: 4 Oct 2024

Richard Camuccio<br>
rcamuccio@gmail.com

(Imageredux > CAL > Hades)