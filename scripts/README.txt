quicktrim.py is a script to subtract bias and trim overscan from images taken with STA1600.

you can display some output for debugging puposes with argument -v

	python quicktrim.py -v

quicktrim needs a bias frame at the project root directory and directories data and trimmed.
The project layout should look as follows:
	\data
	\trimmed
	\bias.fits
	\quicktrim.py
	\README.txt


NOTE: The Bias frame should match the shape of the data frames.
The trimmed frames will be a smaller shape (10560, 10560)

Basic usage
1. drop all files to be bias subtracted and trimmed into the data directory.
2. run quicktrim bash script: right click quicktrim.sh, run as program
	 -or- python directly: python quicktrim.py 
3. get files from trimmed directory

Bash Script
quicktrim.sh exists for ease of use either directly in terminal or
    right click and run as program.

First set quicktrim.sh executable permissions with terminal by running

    chmod +x quicktrim.sh

You can then run quicktrim.sh as a program after setting permissions using
    the instructions in Basic usage.
