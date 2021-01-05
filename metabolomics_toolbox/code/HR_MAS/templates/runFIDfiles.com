#!/bin/csh

# Run the Scripts
	cd $1
	
	foreach filename (`ls`)
        $filename
	end
