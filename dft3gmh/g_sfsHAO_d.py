#!/usr/bin/env python

import numpy as np , os , sys

from scipy import linalg


usage = os.path.basename(__file__) + ": [ symmetric FAO Name ] [ symmetric SAO Name ] [ (optional) HAO Name ]"

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



fao = np.loadtxt( faoName )
sao = np.loadtxt( saoName )
hao = linalg.inv( linalg.sqrtm( sao ) ).dot( fao ).dot( linalg.inv( linalg.sqrtm( sao ) ) )

np.savetxt( haoName , hao ) ;

