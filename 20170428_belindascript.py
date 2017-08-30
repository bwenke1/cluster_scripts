#!/usr/bin/env python
import os 
from sys import argv as argvs 
 
post1 = 1303
post2 = 1456
if post1 > post2:
    print "postion_1 must lower than postion_2, or you will get nothing! Please check it! "

print "Step 1: We will test the files you wanted"
print argvs 

del argvs[0]
if len(argvs) == 0:
    print "No file was pointed out, Please check you input!" 

for files in argvs:
    lines_for_del=[]
    
    if os.path.isfile(files):
        print "================================================================"
        print "The file of:", files, " will be analysised." 
        onefile = open(files)
	for numb, value in enumerate(onefile):
#ok so now it will read the file you put in line by line
	    arrline = value.split() 
# so now each line is split by the spaces.  
	    if len(arrline) > 4: 
#that command skips down to the part where it's not the header
	        post0 = float( arrline[0] )
# now it will convert the value in the 0th position in the split line into a float so you can do the comparison of value
# if you say arrline[0] it is rlnCoordinateX
# if you say arrline[1] it is rlnCoordinateY
	        if post0 >= post1 and post0 <= post2:
	            lines_for_del.append(numb)
	onefile.close()
    else:
       print files,"  is not exist, please check it !"  
       continue 
    
    if len(lines_for_del) > 0:
        oldname = files + '_old'
	os.rename( files, oldname ) 
	oldfile = open(oldname,"r")
        newfile = open(files, "w")
	nline=0
	for numb, value in enumerate(oldfile):
	   if numb not in lines_for_del:
	       newfile.write(value)
	   else:
	       nline+=1
	       print "The number: ", nline, " line will be delected.\n",value
        oldfile.close()
        newfile.close()
        print "The file of:", files, " was renewed" 
	print "****************************************************************"

