#!/bin/csh
# Generates and Runs an ft script (in a simple way)

# Generate the script:

	basicFT1.com -in $1 \
				 -out $3/$4.ft \
				 -xELB 1.5 \
				 -xP0 0 \
				 -xP1 0 \
				 -list > $2/$4_ft.com
	chmod a+rwx $2/$4_ft.com

# More options (order matters):
	# -xELB 1.5 -xP0 146.3
# chmod a+rx ./*.com
	# -in      Read in the fid file
	# -fn SP   Window, (watch first and last points)
	# -fn ZF   Zero Fill (Double Size, Round Up to Power of 2)
	# -fn FT   Fourier Transform
	# -fn PS   Phase Correction
	# -fn TP   Transpose
	# -fn LP   Linear Prediction
	# -fn EXT  Extract Left Half
	# -fn MED  Median Baseline Filter
	# -fn POLY Polynomial Baseline Filter


# Run the Script

	$2/$4_ft.com

# Do referencing (if you want)
#ref1D.tcl -in test.ft1 -x1 0.1ppm -xn -0.1ppm -ref 0.0

# change ft header
#sethdr test.ft1 -title spec
