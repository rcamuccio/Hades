# Hades

**HADES** is an astronomical event listener and data analysis pipeline.

## Installation

### Development platform

Hades was developed on an Alienware Aurora R7 (`Epimetheus`) running Ubuntu 22.04.3 LTS. The general specifications of Epimetheus include an Intel Core i7-8700 x 12 processor, 64 GB of memory, and a disk capacity of 2.3 TB.

Since Hades is written primarily in Python, it is highly recommended to run Hades in a controlled Conda environment. One can load the environment packages using the `hades.yml` file to initialize your environment.

The difference imaging portion of the image analysis pipeline is dependent on an algorithm written in C. One must have the cfitsio package installed on their system. For Ubuntu users, simply run

```
$ sudo apt install libcfitsio-dev
```

### GCN subscription

One requires an account on the NASA General Coordinates Network (GCN) platform (https://gcn.nasa.gov/). The parameters `client_id` and `client_secret` are uniquely generated per user, and must be plugged into the configuration file in order to use the GCN listener.

### Catalogs

The GCN listener uses the GLADE catalog (both the latest GLADE+ and previous Glade 2.4 catalogs, https://glade.elte.hu/) for selecting galaxy targets to follow up. Download the catalogs as ASCII text files (e.g. GLADE+.txt, GLADE_2.4.txt) and save them in a separate directory.

### Astrometry.net

The reduction sequence includes plate solving FITS files. The plate solving routine is handled by a local instance of Astrometry.net (https://astrometry.net/use.html). One needs to install Astrometry.net and the index files it uses for solving fields.

There is a shell script available (https://github.com/rcamuccio/package-installers/blob/master/astrometry-installer.sh/) which can be used to download all of the required index files.

## Usage

The code is run from the main HADES directory through a series of available scripts.

### Configuration

The settings are controlled in the `config.py` file.

The pipeline assumes the following directory format for all data files:

```
[data_directory]
	[bias]
		[yyyy-mm-dd]
			bias.fits
			...
		...
	[calibration]
		[tmp_bias]
		[tmp_dark]
		[tmp_flat]
		bias.fits
		dark.fits
		flat.fits
	[clean]
		[yyyy-mm-dd]
			[FIELD_xx.xxx]
				img.fits
				...
			...
		...
	[darks]
        [yyyy-mm-dd]
            dark.fits
            ...
        ...
	[diff]
        [yyyy-mm-dd]
            [FIELD_xx.xxx]
                img.fits
                ...
            ...
        ...
	[flats]
        [yyyy-mm-dd]
            flat.fits
            ...
        ...
	[flux]
        [yyyy-mm-dd]
            [FIELD_xx.xxx]
                file1.flux
                ...
            ...
	[lc]
        [FIELD_xx.xxx]
            file1.lc
            ...
        ...
	[master]
        [FIELD_xx.xxx]
            [centroids]
            [tmp_master]
                xx_tmp_master.fits
                ...
            FIELD_xx.xxx_master.fits
            FIELD_xx.xxx_star_list.txt
        ...
	[raw]
        [yyyy-mm-dd]
            [FIELD_xx.xxx]
                img.fits
                ...
            ...
        ...
	[review]
        [yyyy-mm-dd]
            [FIELD_xx.xxx]
                img.fits
                ...
            ...
        ...
```

### Running the pipeline

```
$ python -m scripts.main
```

## References

- G. Dálya et al. (2018), GLADE: A Galaxy Catalogue for Multimessenger Searches in the Advanced Gravitational-Wave Detector Era, MNRAS, Volume 479, Issue 2, https://doi.org/10.1093/mnras/sty1703

- G. Dálya et al. (2022), GLADE+: An Extended Galaxy Catalogue for Multimessenger Searches with Advanced Gravitational-Wave Detectors, MNRAS, Volume 514, Issue 1, https://doi.org/10.1093/mnras/stac1443

- LVK Collaboration (2022), LIGO/Virgo/KAGRA Public Alerts User Guide, https://emfollow.docs.ligo.org/userguide/

---

This repository is a development fork of the TOROS pipeline (https://github.com/ryanoelkers/toros/).

June 7, 2019<br>
Last update: January 1, 2026

Richard Camuccio<br>
rcamuccio@gmail.com