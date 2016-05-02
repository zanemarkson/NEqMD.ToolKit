#!/usr/bin/env python

import numpy as np

import os , sys


def gen2g09dftb ( genFile , g09File , paraPath="GAUSS_EXEDIR/dftb/3ob-3-1/" ):

	f = open( genFile , 'r' )

	prefix = g09File[ : -3 ]

	chkFile = prefix + "chk" 

	g = open( g09File , 'w+' )

	lines = f.readlines()

	natom = int( lines[0].split()[0])

	atomTypes = lines[1].split()

	natomTypes = len( atomTypes )

	rawGeom = lines[ 2 : ]

	geom = np.zeros([ natom , 4 ])

	for iatom in np.arange( natom ):

		geom[ iatom , 1 : ] = np.array(rawGeom[ iatom ].split()[ 2 : ] , dtype='float')

		geom[ iatom , 0 ] = rawGeom[ iatom ].split()[1]

	
	g.write("%%chk=%s\n%%mem=4GB\n%%nprocshared=4\n\n#P DFTB Pop=Full IOp(6/7=3,6/80=1) gfinput NoSymm\n\nDFTB+ Parameter\n\n0\t1\n" % chkFile )

	for iatom in np.arange( natom ):

		g.write("  %s\t% 12.8f\t% 12.8f\t% 12.8f\n" % ( atomTypes[ int(rawGeom[ iatom ].split()[ 1 ]) - 1 ] , geom[ iatom , 1 ] , geom[ iatom , 2 ] , geom[ iatom , 3 ] ) )

	g.write("\n")


	for it1 in range( natomTypes ):
		for it2 in range( natomTypes ):
			g.write( "@%s%s-%s.skf/n\n" % ( paraPath , atomTypes[it1] , atomTypes[it2] ) )

	g.write( "\n\n\n" )



usage = os.path.basename(__file__) + ": [DFTB+ File Name] [G09 File Name (to be generated)] [parameter path for G09]"

argc = len( sys.argv )

if argc == 1:
	print "\n%s\n\n" % usage

elif argc == 2 and sys.argv[ 1 ] == "-h":
	print "\n%s\n\n" % usage

elif argc == 4 and sys.argv[ 1 ] != "-h":
	
	genFile = sys.argv[ 1 ]

	g09File = sys.argv[ 2 ]

	paraPath = sys.argv[ 3 ]

	if paraPath[ -1 ] != '/':
		paraPath = paraPath + '/'


	gen2g09dftb ( genFile ,  g09File , paraPath )

elif argc == 3 and sys.argv[ 1 ] != "-h":

	genFile = sys.argv[ 1 ]

	g09File = sys.argv[ 2 ]

	gen2g09dftb ( genFile ,  g09File )

else:

	print("\nExecution Error! Please refere to \"%s -h\"" % os.path.basename(__file__) )



