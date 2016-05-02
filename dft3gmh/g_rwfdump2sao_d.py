#!/usr/bin/env python

import numpy as np , os , sys

from scipy import linalg

def dump2mat( dumpName ):

	with open( dumpName , 'r' ) as f :

		dumpLines = f.readlines()

		for iline in range(len(dumpLines)):

			if (dumpLines[iline].find('read left to right')) == -1:

				continue

			else:

				break

		infoLine = dumpLines[iline]

		nTriElement = int( infoLine.split()[ infoLine.split().index('length') + 1 ] )

		dim = int(np.floor( np.sqrt( nTriElement * 2 ) ))

		rawMatrix = np.array( ' '.join( dumpLines[ iline + 1 : ] ).replace( 'D' , 'E' ).split() ).astype( np.float )

		print("\nrawMatrix %s Dimension is %d\n" % ( dumpName , rawMatrix.shape[0] ) )

		fullMatrix = np.zeros( [ dim , dim ] )

		print("\nDimension of This Matrix is %d-by-%d\n" % ( dim , dim ) )

		id = 0

		for irow in range( dim ):

			for icol in range( irow + 1 ):

				fullMatrix[ irow , icol ] = rawMatrix[ id ]

				fullMatrix[ icol , irow ] = fullMatrix[ irow , icol ]

				id = id + 1

	return fullMatrix


# main function

usage = os.path.basename(__file__) + ": [ RWF-514R Name ] [ (optional) HAO Name ]"

argc = len( sys.argv )

if argc == 1:
	print "\n%s\n\n" % usage
	quit()

elif argc == 2 and sys.argv[ 1 ] == "-h":
	print "\n%s\n\n" % usage
	quit()

elif argc == 2 and sys.argv[ 1 ] != "-h":
	saoName = sys.argv[ 1 ]
	haoName = saoName[ : -3 ] + 'FAO'

elif argc == 3 and sys.argv[ 1 ] != "-h":
	saoName = sys.argv[ 1 ]
	haoName = sys.argv[ 2 ]

else:
	print("\nExecution Error! Please refere to \"%s -h\"" % os.path.basename(__file__) )
	sys.exit( -13 )



sao = dump2mat( saoName )

np.savetxt( haoName , sao ) ;

