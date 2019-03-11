#!/bin/csh

# Arguments:
# 	$1 : relative filepath to directory containing raw files
# 	$2 : relative filepath to directory for processed .fid files
# 	$3 : relative filepath to directory for processed .ft files


# Process the raw .fid files and copy them to a new place, then do nmrPipe on them

	# Navigate to the raw directory (containing one data type), so we can get the file numbers:
			
		cd $1
		
	# Loop through all files in that directory:
		foreach i (`ls`)
		
			# Do initial processing, store resulting .fid:
				bruk2pipe -in $i/fid \bruk2pipeParamsGoHere					-out $2/$i.fid -verb -ov
		end
		
		foreach i (`ls`)		  
			# Make the ft files from the processed .fid files using nmrPipe
				nmrPipe -in $2/$i.fid \
					#| nmrPipe -fn SOL \
					| nmrPipe -fn EM -lb 2 \
					| nmrPipe -fn ZF -auto \
					| nmrPipe -fn FT -auto -verb \
					| nmrPipe -fn PS -p0 0 -p1 -27360 \
					| nmrPipe -fn EXT -x1 -0.5ppm -xn 12.0ppm -sw \
					| nmrPipe -fn PS -p0 0.0 -p1 0.0 -di \
					| nmrPipe -fn POLY -auto\
					| nmrPipe -out $3/$i.ft -ov -verb 
		end

# Remove any empty ft files
	cd $3
	find . -type f -empty -delete
	cd ..







