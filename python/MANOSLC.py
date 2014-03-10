##############################################################
## MANOSLC.py
##
## PURPOSE: Produces .txt file of light curve data in MANOS format
##          and a light curve file in .png format
##
## INPUT: from Polishook iraf script, filenames for refStarsPhot_Final, astMagPhot_Final
##
## OUTPUT: .phot.eps, .manos.txt
##
## HOW TO RUN: called by Polishook iraf script,
##              or 'python MANOSDATA.py file_refStarsPhot_Final file_astPhot_Final from terminal
##
## ANYTHING ELSE (NOTES:
##
## CHANGELOG:
##     2014 01 25 // tge -- combo mask added
##############################################################

import sys
import time
import numpy as np
import matplotlib.pyplot as plt
import os
import os.path

current_dir = os.getcwd() # find the current directory
results_directory ='/Users/manos/catalog/' #permanent results path on allagash
#results_directory = '/Users/a20261/Ureka/python/mytest/manosfiles/' #local dir for debug/testing

#accept command line arguments from Polishook iraf script.
refStarsPhot_Final_file = sys.argv[1] #ex: 2013BO76_20130901_ANDICAM_refStarsPhot_Final
astPhot_Final_file = sys.argv[2] #ex:2013BO76_20130901_ANDICAM_astPhot_Final

#parse filename for obj,instr,date
filename_data = astPhot_Final_file.split('_')
object = filename_data[0]
gregorian_date = filename_data[1]
instrument = filename_data[2]

refStarsPhot_Final_path = current_dir+'/'+refStarsPhot_Final_file #/Users/allagash/Desktop/manosfiles/2013BO76_20130901_ANDICAM_refStarsPhot_Final
astPhot_Final_path = current_dir+'/'+astPhot_Final_file #/Users/allagash/Desktop/manosfiles/2013BO76_20130901_ANDICAM_astPhot_FinalWITHTWOS

#verify object name with user
print 'Is object name "%s" correct? (Y/N)' % (object)
objectok = raw_input()
if objectok.lower() in {'n','no'}:
    object = raw_input('Please enter object designation: ')
#UPDATE NEEDED: add statement to evaluate that rawinput matches MANOS preferred format
elif objectok.lower() in {'y','ye','yes'}:
    pass
else:
    print 'Unexpected object designation, original object name, "%s" will be used.' % (object)



# ASTEROID
#read in astPhot_Final as astobslist and manos_data
with open(astPhot_Final_path,'r') as g:
    astobslist = np.genfromtxt(g,
                               dtype = ["|S30", int, "|S10", float, "|S14", float, float, float],
                               names = ['ImageName','StarNumber','Filter','JulianDate','Airmass','Aperture','Mag','Merr'],
                               usecols = (0,1,3,4,6,7),
                               autostrip=True
                               )
with open(astPhot_Final_path,'r') as g:
    manos_data = np.genfromtxt(g,
                               dtype = ["|S30", int, "|S10", float, "|S14", float, float, float],
                               names = ['ImageName','StarNumber','Filter','JulianDate','Airmass','Aperture','Mag','Merr'],
                               autostrip = True
                               )

#remove duplicate lines from astobslist and manos data
astobslist_nodupes = []
manos_data_nodupes = []
for i in range(0,len(astobslist)):
    if astobslist[i][1]==1:
        astobslist_nodupes.append(astobslist[i])
        manos_data_nodupes.append(manos_data[i])

astobslist = astobslist_nodupes
manos_data = manos_data_nodupes


#create 1D np.arrays for astmags and astmerrs
#with rows equal to # of images
astmags = np.zeros([len(astobslist)])
astmerrs = np.zeros([len(astobslist)])


#populate with corresponding mags/merrs from astobslist
for img_number in range(0,len(astobslist)):
    astmags[img_number] = astobslist[img_number][4]
    astmerrs[img_number] = astobslist[img_number][5]



# REF STARS
#read in refStarsPhot_Final as refobslist
f = open(refStarsPhot_Final_path,'r')
refobslist = np.genfromtxt(f,
                           dtype = ["|S30", int, "|S10", float, "|S14", float, float, float],
                           names = ['ImageName','StarNumber','Filter','Airmass','JulianDate','Aperture','Mag','Merr'],
                           usecols = (0,1,3,4,6,7),
                           autostrip=True
                           )


#create list of img names from asteroid mag file
asteroid_imagenames = tuple(row[0] for row in astobslist)

#retrieve ref mags with matching image names
refobsimgmatches = []
for i in range(0,len(refobslist)):
    if refobslist[i][0] in asteroid_imagenames:
        refobsimgmatches.append(refobslist[i])

#exclude refstar images that do not appear in asteroid image list
refobslist = refobsimgmatches


#calculate total stars in refobslist by finding max value of second column
totalstars=0
for i in range(len(refobslist)):
    if totalstars < refobslist[i][1]:
        totalstars = refobslist[i][1]


num_of_images = len(refobslist)/totalstars
#initialize np.arrays for refmags, refmerrs, and refmeanmags
#with rows equal to # of images, cols equal to # of stars
refmags = np.zeros(shape=(len(refobslist)/totalstars,totalstars))
refmerrs = np.zeros(shape=(len(refobslist)/totalstars,totalstars))


#populate refmags with corresponding magnitudes from refobslist
for img_number in range(0,num_of_images):
    for star in range(0,totalstars):
        refmags[img_number,star] = refobslist[totalstars*img_number+star][4]
        refmerrs[img_number,star] = refobslist[totalstars*img_number+star][5]

#build mask for INDEF magnitude values
#indef_mask = np.ma.MaskedArray(refmags,np.isnan(refmags))
indef_mask = np.zeros(np.shape(refmags))
for img_number in range(0,num_of_images):
    for star in range(0,totalstars):
        if np.isnan(refmags[img_number][star]):
            indef_mask[img_number][star] = 1


#mask for stars with large errors compared to asteroid image
large_err_mask = np.zeros(np.shape(refmags))
for img_number in range(0,num_of_images):
    for star in range(0,totalstars):
        if 3*astmerrs[img_number]<refmerrs[img_number,star]:
            large_err_mask[img_number,star] = 1


#combine indef_mask and large_err_mask
combo_mask = np.zeros(np.shape(refmags))
for img_number in range(0,num_of_images):
    for star in range(0,totalstars):
        if indef_mask[img_number][star] !=0 or large_err_mask[img_number][star] != 0:
                combo_mask[img_number][star] = 1


print 'Prior to median rejection %s of %s stars are good.' \
    % (int(totalstars-np.sum(np.prod(combo_mask,axis=0))), totalstars),'\n'


#create median_mask for stars "far away" from median magnitude
median_mask = np.zeros(np.shape(refmags))
ma_refmags = np.ma.array(refmags,mask=combo_mask)
median_refmags = np.ma.median(ma_refmags,axis=1)
std_refmags = ma_refmags.std(axis=1)

#epsilon gives +/- range of acceptable median values, based on number of stars
#remaining after indef and large_err rejection

# * * * * * * * *
# Review math for epsilon choice as percent of totalstars
# * * * * * * * *
epsilon = totalstars/int(totalstars-np.sum(np.prod(combo_mask,axis=0)))
rejected_stars = [0]*num_of_images #list shows total rejected stars in each img

for img_number in range(0,num_of_images):
    rejected_by_median = 0
    for star in range(0,totalstars):
        if (refmags[img_number][star] <= median_refmags[img_number]-epsilon
            or refmags[img_number][star] >= median_refmags[img_number]+epsilon):
            
            median_mask[img_number][star] = 1
            rejected_by_median += 1
    rejected_stars[img_number]=rejected_by_median
    if rejected_by_median == totalstars:
        print 'All stars in image',astobslist[img_number],'were rejected by median comparison'

#combine combo_mask and median_mask
for img_number in range(0,num_of_images):
    for star in range(0,totalstars):
        if combo_mask[img_number][star] !=0 or median_mask[img_number][star] != 0:
            combo_mask[img_number][star] = 1

print 'After median rejection %s of %s stars are good.' \
    % (int(totalstars-np.sum(np.prod(combo_mask,axis=0))), totalstars),'\n'

#create masked refmags and refmerrs arrays
ma_refmags = np.ma.array(refmags,mask=combo_mask)
ma_refmerrs = np.ma.array(refmerrs,mask=combo_mask)


best_image_number = np.argmin(ma_refmerrs.mean(axis=1))
print 'Image used as reference image: ', \
    refobslist[totalstars*best_image_number][0]

best_image_merrs = ma_refmerrs[best_image_number]
ma_refmags_diff = ma_refmags - ma_refmags[best_image_number]
ext_corr_astmags = astmags - ma_refmags_diff.mean(axis=1)

#calculate final asteroid magnitude, centered at zero
final_astmags = np.around(ext_corr_astmags - ext_corr_astmags.mean(),5)

#calcualte final asteroid magnitude error

# * * * * * * * * * * *
# -- possible error calculations --
# 1.
final_astmerrs = np.around(((astmerrs**2)+ \
                            (ma_refmags_diff.std(axis=1))**2+ \
                            (ma_refmerrs.mean(axis=1))**2)**.5,5)
# 2.
#final_astmerrs = np.around(((astmerrs**2)+(ma_refmerrs.mean(axis=1))**2)**.5,5)

# 3.
#final_astmerrs = np.around(((astmerrs**2)+
#                            (ma_refmerrs.mean(axis=1))**2+
#                            ((ma_refmerrs - best_image_merrs).mean(axis=1))**2)**.5,5)
# 4.
#final_astmerrs = np.around(((astmerrs**2)+(ma_refmags_diff.std(axis=1))**2+((ma_refmerrs - best_image_merrs).mean(axis=1))**2)**.5,5)

# * * * * * * * * * * *

#update manos_data array with calibrated magnitudes
for i in range(0,len(manos_data)):
    manos_data[i][6] = final_astmags[i]
    manos_data[i][7] = final_astmerrs[i]

## remove lines with INDEF in mag/merr
#manos_data[~np.isnan(final_astmags).any(axis=1)]
#manos_data[~np.isnan(final_astmerrs).any(axis=1)]
##


#check if directory results_directory/OBJECT exists, if not, create directory
if not os.path.isdir(results_directory+'n_'+object+'/') \
        and not os.path.isdir(results_directory+object+'/'):
    if object.isdigit():
        results_directory = results_directory+'n_'+object+'/'
    else:
        results_directory = results_directory+object+'/'
    os.makedirs(results_directory)
#if directory does exist, update the results_directory variable with object name
else:
    if object.isdigit(): results_directory = results_directory+'n_'+object+'/'
    else: results_directory = results_directory+object+'/'


#check if object is a numbered object, add 'n_' to filename
if object.isdigit():
    fileout = 'n_'+object+'_'+gregorian_date
else:
    fileout = object+'_'+gregorian_date

#check if manos.txt, phot.eps exist, if so add next sequence number
segment = 1
while os.path.isfile(results_directory+fileout+'.manos.txt') \
    or os.path.isfile(results_directory+fileout+'.'+str(segment)+'.manos.txt'):
    segment += 1
fileout = fileout+'.'+str(segment)


#save manos_data as .manos.txt
with open(results_directory+fileout+'.manos.txt','w') as h:
    np.savetxt(h,manos_data,delimiter=' ', fmt='%6s')

#PLOTTING ROUTINE

#declare three empty lists to hold date, mag, merrs
#find middle JD in astobs for MPC data
start_JD = '%.4f' % manos_data[0][3]
mid_JD = str(int(manos_data[len(manos_data)/2][3]))
x_hours=[]
y_mags=[]
y_err=[]

for i in range(0,len(manos_data)):
    if str(manos_data[2]) != 'INDEF' or np.isnan(manos_data[2]):
        x_hours.append(25*(-.5+manos_data[i][3]%1))
        y_mags.append(manos_data[i][6])
        y_err.append(manos_data[i][7])

minimum_magnitude = np.argmin(y_mags)
maximum_hour = np.argmax(x_hours)
minimum_hour = np.argmin(x_hours)

#display plot and save as .phot.eps
plt.plot(x_hours-x_hours[0], y_mags, 'ro')
plt.errorbar(x_hours-x_hours[0], y_mags, yerr = y_err, ls='none', color='black')
plt.tick_params(axis='both',labelsize=10)
plt.xlim([min(x_hours)-x_hours[0]-.05*(max(x_hours)-min(x_hours)),
          max(x_hours)-x_hours[0]+.05*(max(x_hours)-min(x_hours))])
plt.xlabel('Hour')
plt.ylabel('Magnitude')
plt.title(object+'   '+instrument+'   '+gregorian_date)
plt.gca().invert_yaxis()
plt.annotate('Start Time:\n'+start_JD,ha='left',size=10,
                xy=(0.025,0.025), xycoords='axes fraction')
plt.savefig(results_directory+fileout+'.phot.eps',format='eps')

print 'Light curve and MANOS file saved in \n',results_directory

#open .eps  file (mac os)
os.system('open '+results_directory+fileout+'.phot.eps')

#call MPCDATA script using .manos.txt
#import creatempcdatafile
#creatempcdatafile.writefile(fileout+'.manos.txt',0.0)


#
#
#

