#!/bin/bash
# $1	 relative path from pipe_templates to directory containing fid files
		# we will read in the filenames under this directory, and pass them to sed
# $2	 string to replace (relative path to fid file)
# $3	 relative path/input file that we will operate on and copy to _ft.com using sed
# $4	 relative path to ft.com output destination

# Navigate from this file (template dir) to the fid files directory:
	cd $1
	
# Loop through the filenames in the directory
	ls | while read fname ; do 
	# Test that we get the filenames:
	
		#echo "${fname%.*}";	# cut off the extension
		#echo "in file $3, replace all instances of $2 with ${fname%.*}, then store the output in $4/${fname%.*}_ft.com"
		#echo "${fname%.*}"
		
	# Test the sed substitution command to replace the existing filename
	# with the current fname:
	
		#echo "'s/$2/${fname%.*}/g' $3 > $4/${fname%.*}_ft.com;"
		
	# Run the command
########	
		#specNumber = "${fname%.*}" # add later for speed
########	

		sed "s/$2/${fname%.*}/g" $3 > $4/${fname%.*}_ft.com;
		
	# Make the file executable
	
		chmod a+x $4/${fname%.*}_ft.com;
		
########			
    # Run the file
        
        #echo "$4/${fname%.*}_ft.com"
        #$4/${fname%.*}_ft.com			# add later for speed
########	

	done

