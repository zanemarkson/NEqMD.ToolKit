#!/usr/bin/env python

import numpy as np, os, sys

from scipy import linalg


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


def loadEnergyFile(energyFileName):
    with open(energyFileName, 'r') as f:
        dumpLines = f.readlines()
        MOEnergy = np.array(' '.join(dumpLines).replace('D', 'E').split()).astype(np.float)
        return MOEnergy


def readpyci(cifilename):
    with open(cifilename, 'r') as f:
        allLines = f.readlines()

        f.seek(0)
        gsLine = allLines[fsearch(f, 'GS')].split()[2:]
        gs = np.array(gsLine).astype(np.int)
        # n_gs = gs.shape[0]

        f.seek(0)
        ctLine = allLines[fsearch(f, 'CT')].split()[2:]
        ct = np.array(ctLine).astype(np.int)
        # n_ct = ct.shape[0]

        f.seek(0)
        leLine = allLines[fsearch(f, 'LE')].split()[2:]
        le = np.array(leLine).astype(np.int)
        # n_le = le.shape[0]

        return gs, ct, le


def whichMO(MOs, candidateList, aoList, nselect):
    nCandidate = len(candidateList)
    performance = np.zeros(nCandidate, dtype={'names': ['id', 'perform'], 'formats': ['i4', 'f4']})
    performance['id'] = candidateList
    for iCandidate in np.arange(nCandidate):
        thisCandidate = candidateList[iCandidate]
        thisPerformance = np.sum(MOs[aoList, thisCandidate]**2)
        performance[iCandidate]['perform'] = thisPerformance
    performance.sort(order='perform')
    return performance[::-1]['id'][0:nselect]


def gmh(approx_type, dipX, dipY, dipZ, energies, MO, gsOrbs, ctOrbs, leOrbs):
    # AO Dipole -> MO Dipole

    MODipX = (MO.T).dot(dipX).dot(MO)
    MODipY = (MO.T).dot(dipY).dot(MO)
    MODipZ = (MO.T).dot(dipZ).dot(MO)

    # GS-CT First:
    igsct = 0
    gsct = np.zeros((len(gsOrbs) * len(ctOrbs), 7 ))
    for gs in gsOrbs:
        for ct in ctOrbs:
            h11 = energies[gs - 1]
            h22 = energies[ct - 1]
            dE = np.abs(h11 - h22)
            mu11 = np.array([MODipX[gs - 1, gs - 1], MODipY[gs - 1, gs - 1], MODipZ[gs - 1, gs - 1]])
            mu22 = np.array([MODipX[ct - 1, ct - 1], MODipY[ct - 1, ct - 1], MODipZ[ct - 1, ct - 1]])
            mu12 = np.array([MODipX[gs - 1, ct - 1], MODipY[gs - 1, ct - 1], MODipZ[gs - 1, ct - 1]])
	    print "\n{0}->{1} Transition Dipole is X: {2} , Y: {3} , Z: {4}".format(gs,ct,MODipX[gs - 1, ct - 1], MODipY[gs - 1, ct - 1], MODipZ[gs - 1, ct - 1])
            dMu = mu11 - mu22
            v0 = dMu / np.linalg.norm(dMu)
            if approx_type == 1:
                mu12v = mu12.dot(v0)
            elif approx_type == 2:
                mu12v = np.linalg.norm(mu12)
            gsct[igsct, 0] = np.linalg.norm(mu12v) * dE / np.sqrt(
                np.linalg.norm(dMu) ** 2 + 4.00 * np.abs(mu12v) ** 2) * 27.211
            gsct[igsct, 1] = dE
            gsct[igsct, 2] = mu12v
            gsct[igsct, 3:6] = dMu
            gsct[igsct, 6] = np.linalg.norm(dMu)
            igsct += 1

    # LE-CT Second:
    ilect = 0
    lect = np.zeros((len(leOrbs) * len(ctOrbs), 7 ))
    for le in leOrbs:
        for ct in ctOrbs:
            h11 = energies[le - 1]
            h22 = energies[ct - 1]
            dE = np.abs(h11 - h22)
            mu11 = np.array([MODipX[le - 1, le - 1], MODipY[le - 1, le - 1], MODipZ[le - 1, le - 1]])
            mu22 = np.array([MODipX[ct - 1, ct - 1], MODipY[ct - 1, ct - 1], MODipZ[ct - 1, ct - 1]])
            mu12 = np.array([MODipX[le - 1, ct - 1], MODipY[le - 1, ct - 1], MODipZ[le - 1, ct - 1]])
	    print "\n{0}->{1} Transition Dipole is X: {2} , Y: {3} , Z: {4}".format(le,ct,MODipX[le - 1, ct - 1], MODipY[le - 1, ct - 1], MODipZ[le - 1, ct - 1])
            dMu = mu11 - mu22
            v0 = dMu / np.linalg.norm(dMu)
            if approx_type == 1:
                mu12v = mu12.dot(v0)
            elif approx_type == 2:
                mu12v = np.linalg.norm(mu12)
            lect[ilect] = np.linalg.norm(mu12v) * dE / np.sqrt(
                np.linalg.norm(dMu) ** 2 + 4.00 * np.abs(mu12v) ** 2) * 27.211
            lect[ilect, 1] = dE
            lect[ilect, 2] = mu12v
            lect[ilect, 3:6] = dMu
            lect[ilect, 6] = np.linalg.norm(dMu)
            ilect += 1
    return gsct, lect


# main function

usage = os.path.basename(__file__) + ": [ rwf-dipole File Name ] [ rwf-MO File Name ] [ Orbital Energy File ] " \
                                     "[ pyCI file name ][ (optional) coupling file Name ] [ optional dipole approx type]"

argc = len(sys.argv)

if argc == 1:
    print "\n%s\n\n" % usage
    quit()

elif argc == 2 and sys.argv[1] == "-h":
    print "\n%s\n\n" % usage
    quit()

elif (argc == 3 or argc == 4) and sys.argv[1] != "-h":
    print "\n{0}\n".format('Error! Please refer to -h')
    quit()

elif argc == 5 and sys.argv[1] != "-h":
    dipName = sys.argv[1]
    moName = sys.argv[2]
    energyName = sys.argv[3]
    pyCIName = sys.argv[4]
    bdiagName = moName[: -3] + 'bdiag'
    dip_approx_type = 2


elif argc == 6 and sys.argv[1] != "-h":
    dipName = sys.argv[1]
    moName = sys.argv[2]
    energyName = sys.argv[3]
    pyCIName = sys.argv[4]
    bdiagName = sys.argv[5]
    dip_approx_type = 2

elif argc == 7 and sys.argv[1] != "-h":
    dipName = sys.argv[1]
    moName = sys.argv[2]
    energyName = sys.argv[3]
    pyCIName = sys.argv[4]
    bdiagName = sys.argv[5]
    dip_approx_type = int(sys.argv[6])

else:
    print("\nExecution Error! Please refere to \"%s -h\"" % os.path.basename(__file__))
    sys.exit(-13)

homo = 390 - 1
AO_ccAntro = np.arange(1021, 1290) - 1
AO_dmaAndG = np.arange(1, 345) - 1


dipX, dipY, dipZ = dump2dip(dipName)
MO = dump2mo(moName)
energies = loadEnergyFile(energyName)

# ===> Here we added 1 to match the gs, ct, le definition inside function gmh()
gs = whichMO(MO, np.arange(homo, homo-4, -1), AO_ccAntro, 1) + 1
ct = whichMO(MO, np.arange(homo, homo-4, -1), AO_dmaAndG, 2) + 1
le = np.array([390]) + 1

print "\ngs = [ {0} ]\t;\tct = [ {1} ]\t;\tle = [ {2} ]\n".format(gs, ct, le)

# gs, ct, le = readpyci(pyCIName)

gsct, lect = gmh(dip_approx_type, dipX, dipY, dipZ, energies, MO, gs, ct, le)

couplings = np.concatenate((gsct[:,0], lect[:,0]))

dEs = np.concatenate((gsct[:,1], lect[:,1]))

transitionDipoleMomentNorms = np.concatenate((gsct[:,2], lect[:,2]))

differenceDipoleMomentNorms = np.concatenate((gsct[:,6], lect[:,6]))

rendering = np.concatenate((couplings, dEs, transitionDipoleMomentNorms,
                                                          differenceDipoleMomentNorms))

np.savetxt(bdiagName, rendering, fmt='% .10E', delimiter='\t')




MODipX = (MO.T).dot(dipX).dot(MO)
MODipY = (MO.T).dot(dipY).dot(MO)
MODipZ = (MO.T).dot(dipZ).dot(MO)

# np.savetxt('MODipX.deb', MODipX, fmt='% .10E', delimiter='\t' )
# np.savetxt('MODipY.deb', MODipY, fmt='% .10E', delimiter='\t' )
# np.savetxt('MODipZ.deb', MODipZ, fmt='% .10E', delimiter='\t' )
print "\nTrace of MO Dipole Matrices are:\ntr[ MODipX ] = [ {0} ]\t;\n\ntr[ MODipY ] = [ {1} ]\t;\n\ntr[ MODipZ ] = [ {2} ]\n\n".format(np.trace(MODipX), np.trace(MODipY), np.trace(MODipZ))

# np.savetxt('AOMO.deb', MO , fmt='% .10E', delimiter='\t' )
# np.savetxt('MOEnergy.deb', energies , fmt='% .10E', delimiter='\n' )



