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


def tri2sym(rawMatrix):
    nTriElement = len(rawMatrix)
    dim = int(np.floor(np.sqrt(nTriElement * 2)))

    fullMatrix = np.zeros([dim, dim])
    id = 0
    for irow in range(dim):
        for icol in range(irow + 1):
            fullMatrix[irow, icol] = rawMatrix[id]
            fullMatrix[icol, irow] = fullMatrix[irow, icol]
            id += 1
    return fullMatrix



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



def dump2dip(dumpName):
    with open(dumpName, 'r') as f:
        dumpLines = f.readlines()
        f.seek(0)
        iline = fsearch(f, 'read left to right')

        infoLine = dumpLines[iline]
        nAllElement = int(infoLine.split()[infoLine.split().index('length') + 1])
        nTriElement = nAllElement / 3

        # dim = int(np.floor(np.sqrt(nTriElement * 2)))
        # print("\nrawMatrix %s Dimension is %d\n" % (dumpName, rawMatrix.shape[0]))
        rawMatrixAll = np.array(' '.join(dumpLines[iline + 1:]).replace('D', 'E').split()).astype(np.float)

        dipX = tri2sym(rawMatrixAll[0 * nTriElement:1 * nTriElement])
        dipY = tri2sym(rawMatrixAll[1 * nTriElement:2 * nTriElement])
        dipZ = tri2sym(rawMatrixAll[2 * nTriElement:3 * nTriElement])

        return dipX, dipY, dipZ






# main function

usage = os.path.basename(__file__) + ": [ RWF-524R Name ] [ (optional) Dipole File Prefix ]"

argc = len( sys.argv )

if argc == 1:
	print "\n%s\n\n" % usage
	quit()

elif argc == 2 and sys.argv[ 1 ] == "-h":
	print "\n%s\n\n" % usage
	quit()

elif argc == 2 and sys.argv[ 1 ] != "-h":
	rwfName = sys.argv[ 1 ]
	dipPrefix = 'None'
	aoDipXName = "AODipX.deb"
	aoDipYName = "AODipY.deb"
	aoDipZName = "AODipZ.deb"

elif argc == 3 and sys.argv[ 1 ] != "-h":
	rwfName = sys.argv[ 1 ]
	dipPrefix = sys.argv[ 2 ]
	aoDipXName = dipPrefix + ".dipX"
	aoDipYName = dipPrefix + ".dipY"
	aoDipZName = dipPrefix + ".dipZ"

else:
	print("\nExecution Error! Please refere to \"%s -h\"" % os.path.basename(__file__) )
	sys.exit( -13 )



ao_dip_X , ao_dip_Y , ao_dip_Z = dump2dip( rwfName )

np.savetxt( aoDipXName, ao_dip_X, fmt='% .10E', delimiter='\t' )
np.savetxt( aoDipYName, ao_dip_Y, fmt='% .10E', delimiter='\t' )
np.savetxt( aoDipZName, ao_dip_Z, fmt='% .10E', delimiter='\t' )


