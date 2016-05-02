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

usage = os.path.basename(__file__) + ": [ RWF-536R Name ] [ RWF-514R Name ] [ (optional) HAO Name ]"

argc = len( sys.argv )

if argc == 1:
	print "\n%s\n\n" % usage
	quit()

elif argc == 2 and sys.argv[ 1 ] == "-h":
	print "\n%s\n\n" % usage
	quit()

elif argc == 3 and sys.argv[ 1 ] != "-h":
	faoName = sys.argv[ 1 ]
	saoName = sys.argv[ 2 ]
	haoName = faoName[ : -3 ] + 'HAO'

elif argc == 4 and sys.argv[ 1 ] != "-h":
	faoName = sys.argv[ 1 ]
	saoName = sys.argv[ 2 ]
	haoName = sys.argv[ 3 ]

else:
	print("\nExecution Error! Please refere to \"%s -h\"" % os.path.basename(__file__) )
	sys.exit( -13 )



fao = dump2mat( faoName )
sao = dump2mat( saoName )
hao = linalg.inv( linalg.sqrtm( sao ) ).dot( fao ).dot( linalg.inv( linalg.sqrtm( sao ) ) )

#h = open( haoName , 'w+' )

nrow = hao.shape[0]
ncol = hao.shape[1]

np.savetxt( haoName , hao ) ;
# np.savetxt( 'fao.deb' , fao ) ;
# np.savetxt( 'sao.deb' , sao ) ;

#for irow in range( nrow ):
#	for icol in range( ncol ):
#		h.write( "% 16.12E\t" % hao[ irow , icol ] )
#	h.write( "\n" )
