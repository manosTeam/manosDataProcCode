##############################################################
## creatempcdatafile.py
##
## PURPOSE: Produces a .txt file of light curve data
## in MPC format from manos formatted data
##
## INPUT: .manos.txt file, offset(optional) (eg 'OBJECT_DATE.manos.txt 3.2')
## OUTPUT: .mpc.txt
##
## HOW TO RUN: called from command line
##       as 'python creatempcdatafile.py file.manos.txt offset'
##
##
## CHANGELOG:
##      2014/01/17 // tge -- removed 'DATA' from header information
##                           added for loop to prepend 'DATA=' to all data lines
##      2014/01/23 // tge -- added filter check for MPC allowed filters
##                        -- url name lookup for numbered objects
##                        -- logic for inclusion/exclusion of mpc_desig,object_name
##      2014/01/27 // tge -- added directory/subdirectory logic to store final result
##############################################################
import sys
import numpy as np
import fileinput
import time
import urllib


manos_file = sys.argv[1]
instrument = sys.argv[2]



def writefile(manos_file, offset=0.0):
    
    observer_dictionary = {
            'ANDICAM':'Moskovitz, Nicholas; Polishook, David',
            'TESTSCOPE':'Endicott, Thomas; Moskovitz, Nicholas',
            'SOAR':'Moskovitz, Nicholas',
            'KITTPEAK':'Moskovitz,Nicholas',
            'AstPhot':'Moskovitz,Nicholas',
            'WiseTest':'Polishook,David'
                            }
    filename_data = manos_file.split('_')
    numbered_object = False
    # Commented out telescope variable, instead use general manos contact info
    #telescope = instrument
    results_directory ='/Users/manos/catalog/'       #results path
    url = 'http://www.minorplanetcenter.net/iau/lists/MPNames.html'  #MPC numbered obj list


    #update results directory with object sub directory
    if filename_data[0] != 'n':
        results_directory = results_directory+filename_data[0]+'/'
    else:
        results_directory = results_directory+filename_data[0]+'_'+filename_data[1]+'/'

    # if numbered object, null result for mpc_desig, and lookup object_name
    # else null result for object_name, and mpc_desig == object_number
    if filename_data[0] == 'n':
        filename_data = filename_data[1:]
        numbered_object = False
        object_number = (filename_data[0])
        for line in urllib.urlopen(url).readlines():
            if '('+object_number+')' in line.split():
                object_name = line.split(') ')[1].split('\n')[0].split('<')[0]
                mpc_desig = "" #null value for mpc_desig
                numbered_object = True
        if numbered_object != True:
            print 'Object number does not appear in MPC numbered object list'
            print 'at',url,'\n','Please enter MPC designation in the form \'YYYY XX##\' \n'
            object_name = ""
            user_provided_mpc_desig = raw_input()
            object_number = user_provided_mpc_desig
            mpc_desig = '\n'+'MPCDESG='+object_number
    else:
        object_number = ""
        object_name = filename_data[0][:4]+' '+filename_data[0][4:]
        mpc_desig = '\n'+'MPCDESG='+object_name

    #observers = observer_dictionary[telescope]
    #ignore 'observers' dictionary, used general fixed value instead
    observers = 'The Mission Accessible Near-Earth Asteroid Survey (MANOS)'

    with open(results_directory+manos_file,'r') as manos:
                            manos_data = np.genfromtxt(results_directory+manos_file,
                            dtype = ["|S10", int, "|S10", float, float, float, float, float],
                            names = ['ImageName','StarNumber','Filter','JulianDate','Airmass','Aperture','Mag','Merr'],
                            usecols = (3,6,7,4,2),
                            autostrip = True
                            )
    #convert Julian Date to Gregorian Date
    JD = manos_data[len(manos_data)/2][0]
    if JD>0:
        JDint = int(JD+.5)
    else:
        JDint = int (JD-.5)
    W = int((JDint - 1867216.25)/36524.25)
    X = int(W/4.)
    A = JDint+1+W-X
    B = A+1524
    C = int((B-122.1)/365.25)
    D = int(365.25*C)
    E = int((B-D)/30.6001)
    F = int(30.6001*E)
    day = int(B-D-F+int((JD-JDint)))
    if (E-1)<13:
        month = E-1
    else:
        month = E-13
    if month<3:
        year = int(C-4715)
    else:
        year = int(C-4716)

    #rough check that Gregorian date is a valid date
    #create filename template based on object and date
    if (day>=1 and day<=31) and (month>=1 and month<=12) and year<=time.strftime('%Y'):
        session_date = '%04d-%02d-%02d' % (year, month, day)

    filter_dictionary = {"B":"B","V":"V","R":"R","I":"I","C":"C",
            "u'":"SU","g'":"SG","r'":"SR","i'":"SI","z'":"SZ"
            }

    if manos_data[0][4] in filter_dictionary.values():
        filter = filter_dictionary[manos_data[0][4]]
    else:
        print 'Filter:',manos_data[0][4],
        print 'in MANOS data file is not allowed by MPC guidelines.'
        print 'Allowed filters are:'
        print """B  V  R  I  C  u'  g'  r'  i'  z'"""
        print 'Please provide filter used for observations:'
        filter_provided = raw_input()
        if filter_provided in filter_dictionary:
            filter = filter_dictionary[filter_provided]
        else:
            print 'Provided filter is not allowed by MPC submission guidelines'
            print 'MPC Data file not created. \n End program.'
            sys.exit()

    midJD_decimal = (manos_data[len(manos_data)/2][0]-.5)%1
    session_hour = midJD_decimal*24
    session_minute = (session_hour%1)*60
    session_second = (session_minute%1)*60
    session_time = '%02d:%02d:%02d' % (session_hour,session_minute,session_second)



    template = '''STARTMETADATA
CONTACTINFO=allagash@manos.edu
CONTACTNAME=N. Moskovitz
FILTER={filter}
LTCAPP=NONE
LTCTYPE=None
DELIMITER=PIPE
DIFFERMAGS=TRUE
MAGBAND=V{mpc_desig}
OBJECTNAME={object_name}
OBJECTNUMBER={object_number}
OBSERVERS={observers}
REDUCEDMAGS=NONE
REVISEDDATA=FALSE
SESSIONDATE={session_date}
SESSIONTIME={session_time}
ENDMETADATA
'''
    context = {
        'filter':filter,
        'mpc_desig':mpc_desig,
        'object_name':object_name,
        'object_number':object_number,
        'observers':observers,
        'session_date':session_date,
        'session_time':session_time
        }

    with open(results_directory+manos_file,'r') as f:
        mpc_data = np.genfromtxt(f,
                           dtype = ["|S30", int, "|S10", float, "|S14", float, float, float],
                           names = ['ImageName','StarNumber','Filter','JulianDate','Airmass','Aperture','Mag','Merr'],
                           usecols = (3,6,7,4),
                           autostrip=True
                           )
    for j in range(0,len(mpc_data)):
        mpc_data[j][1] = mpc_data[j][1]+offset

    mpc_file = manos_file[:-10]+'.mpc.txt'

    with open(results_directory+mpc_file,'w') as g:
        g.write(template.format(**context))
    with open(results_directory+mpc_file,'a') as h:
        np.savetxt(h,mpc_data,delimiter='|',fmt='%s')
        h.write('ENDDATA')
    for line in fileinput.input(results_directory+mpc_file, inplace=1):
        if line[0].isdigit():
            print 'DATA='+line,
        else:
            print line,

    print 'MPC data file saved in \n',results_directory


if __name__ == '__main__':
    writefile(manos_file,offset=0.0)

#MPCDESIG is optional, example '1999 CZ1', see OBJECTNAME documentation also
#COMMENT is optional example 'COMMENT=The asteroid was in a very crowded field.'

#FILTER allowed values: B,V,R,I,C(for Clear/No),SU,SG,SR,SI,SZ(case-insensitive)
#   R,I understood as Rc, Ic, indicate if elsewise in COMMENT
#OBJECTNAME name assigned by IAU - if none, enter MPC designation here and for MPCDESIG keyword
#OBJECTNUMBER num assigned by MPC to asteroid, if none value should be empty
#OBSERVERS names of those working telescope and measured images, multiple names separated by ; Lastname, Joe; Doe, Jane
#MPCDESIG={mpc_desig} - optional keyword, not included, prioritize OBJ name and number
#REVISEDDATA - indicates if subimssion should replace previous data
#SESSIONDATE in format of YYYY-MM-DD (approx mid UT date of dataset)
#SESSIONTIME in format HH:MM:SS (approx mid UT time of dataset)
#DATA=JD|MAG|MERR|AIRMASS


#else read in default MANOS metadata from permanent url"

