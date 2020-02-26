#! /usr/bin/env python
# path the to the directory must be provided as command-line argument

import os,sys
for filename in os.listdir(sys.argv[1]):
	if ' ' in filename:
		temp=filename.split(' ')
		temp=''.join(temp)
		os.rename(sys.argv[1]+'/'+filename, sys.argv[1]+'/'+temp)

