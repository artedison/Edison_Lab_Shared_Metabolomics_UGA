#!/bin/csh

# Arguments:
# 	$1 : relative filepath from this file to fid files
# 	$2 : relative filepath from fid files to ft files
#	$3 : relative filepath from this file to ft files

# Process the raw .fid files and copy them to a new place, then do nmrPipe on them

	# Navigate to the raw directory (containing one data type), so we can get the file numbers:
		thisFile = PWD	
		cd $1
				
		foreach i (`ls`)		  
			# Make the ft files from the processed .fid files using nmrPipe
				# Thisfile: pipe_templates
				# File in: from specnum_ft.com in rep spec
				# 
				# In: 	specNum.fid in proc/.../fid 
				# Out: 	specNum.ft in proc/.../ft

				nmrPipe -in ./$i.fid \
					#| nmrPipe -fn SOL \
					| nmrPipe -fn EM -lb 2 \
					| nmrPipe -fn ZF -auto \
					| nmrPipe -fn FT -auto -verb \
					| nmrPipe -fn PS -p0 0 -p1 -0 \
					| nmrPipe -fn EXT -x1 -0.5ppm -xn 12.0ppm -sw \
					| nmrPipe -fn PS -p0 0.0 -p1 0.0 -di \
					| nmrPipe -fn POLY -auto\
					| nmrPipe -out $2/$i.ft -ov -verb 
		end

# Remove any empty ft files
	cd $3
	find . -type f -empty -delete
	cd thisFile


