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


def fsearch(pfile, pattern):
    pfile.seek(0)
    iline = 0
    lines = pfile.readlines()
    for line in lines:
        if (line.find(pattern)) == -1:
            iline += 1
            continue
        else:
            break
    return iline


def dump2mo(dumpName):
    with open(dumpName, 'r') as f:
        dumpLines = f.readlines()
        f.seek(0)
        iline = fsearch(f, 'read left to right')
        infoLine = dumpLines[iline]
        nAllElement = np.int(infoLine.split()[infoLine.split().index('length') + 1])
        dim = np.int(np.sqrt(nAllElement))
        # print("\nrawMatrix %s Dimension is %d\n" % (dumpName, rawMatrix.shape[0]))
        rawMatrixAll = np.array(' '.join(dumpLines[iline + 1:]).replace('D', 'E').split()).astype(np.float)
        MO = rawMatrixAll.reshape(dim, dim)

        return MO.T







# main function

usage = os.path.basename(__file__) + ": [ RWF-536R Name ] [ (optional) MO File Name ]"

argc = len( sys.argv )

if argc == 1:
	print "\n%s\n\n" % usage
	quit()

elif argc == 2 and sys.argv[ 1 ] == "-h":
	print "\n%s\n\n" % usage
	quit()

elif argc == 2 and sys.argv[ 1 ] != "-h":
	rwfName = sys.argv[ 1 ]
	moName = rwfName[ : -3 ] + 'MO'

elif argc == 3 and sys.argv[ 1 ] != "-h":
	rwfName = sys.argv[ 1 ]
	moName = sys.argv[ 2 ]

else:
	print("\nExecution Error! Please refere to \"%s -h\"" % os.path.basename(__file__) )
	sys.exit( -13 )



mo = dump2mo( rwfName )

np.savetxt( moName , mo ) ;

