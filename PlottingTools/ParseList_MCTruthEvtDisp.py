# Python script to parse list of event information and automatically produce MCTruth event displays. Relies very heavily on the exact format of the event information text file, so this code will need to change if that does.

# Takes two arguments (and always needs both). Run as:
#
# python name_of_input_text_file_containing_event_information.txt output_directory_name
#
# If the output directory does not exist, the script will make it.

# Note that this can take a long time for a large number of events! Best to have running in a screen session.

import sys, subprocess
import os, time
from datetime import datetime, timedelta
import signal
import termios
import string
import urllib, urllib2
import ROOT
import glob

inputFileName = sys.argv[1]
outdir = sys.argv[2]

print "Reading events from file " + inputFileName
print "Writing event displays to " + outdir

if not os.path.exists(outdir):
    print "Making output directory..."
    os.mkdir(outdir)

inputFile = open(inputFileName,"r")

input_lines = inputFile.readlines();

outputFileName = outdir+"/evdinfo_filenames.txt"
outputFile = open(outputFileName,"w")

#Loop through lines in file and get event information
events_processed = 0

run = 0
subrun = 0
event = 0
xpos = 0
ypos = 0
zpos = 0
rootfile = "blah"
for iline in range(0,len(input_lines)):
    if events_processed >= 100:
        break
    line = input_lines[iline]
    if ("Run/Subrun/Event:" in line):
        # Set run and event number
        run = line.split()[1]
        subrun = line.split()[2]
        event = line.split()[3]
        # Find root file for this event
        sam_file = subprocess.Popen("samweb list-files \"defname:prodgenie_bnb_nu_cosmic_uboone_mcc8.7_reco2_dev and run_number = %s.%s and first_event <= %s and last_event >= %s\""%(run,subrun,event,event), shell=True, stdout=subprocess.PIPE).stdout.read()
        rootfile = subprocess.Popen("grep %s /pnfs/uboone/persistent/users/ddevitt/filelist/stage1_out_files.txt"%sam_file[:-6], shell=True, stdout=subprocess.PIPE).stdout.read()
        #print rootfile[:-1]
        # Write the rootfile name to the new file so we don't have to look it up
        outputFile.write("------------------\n")
        outputFile.write(rootfile)
        outputFile.write("\n\n")
    if ("Reco nu vertex position (x,y,z):" in line):
        # Set x/y/z positions
        xpos = line.split()[5]
        ypos = line.split()[6]
        zpos = line.split()[7]
        # Now run script to make plot
        rootcommand = "root -x -b -q \"/uboone/app/users/kduffy/CC1pi/CC1pi_uboonecode/srcs/uboonecode/uboone/CC1pi/PlottingTools/mctruth_evtdisplay.C(\\\"%s\\\",%s,%s,%s,%s,\\\"%s\\\")\""%(rootfile[:-1],event,xpos,ypos,zpos,outdir)
        #print rootcommand
        os.system(rootcommand)
        # Finally, reset all values before carrying on in the loop over lines in the file
        run = 0
        event = 0
        xpos = 0
        ypos = 0
        zpos = 0
        events_processed += 1
        rootfile = "blah"
    # Write line to new file
    outputFile.write(line)
