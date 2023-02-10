import numpy as np
import pandas as pd


def find_octave_all_matrices(sourceFilePath):
    KNOWN_OCTAVE_TYPES = {
        # typeName: headerLength
        'scalar': 2,
        'matrix': 4,
        'complex matrix': 4,
    }
    with open(sourceFilePath, 'r') as sourceFile:
        matList = []
        # TODO: Remove following line if Variant 2 is adopted
        # matIndToCheckEnd = -1
        lineInd = 0
        for line in sourceFile:
            lineInd += 1
            if '# name:' in line:
                startInd = lineInd
                typeLine = next(sourceFile)
                if '# ndims:' in next(sourceFile):
                    ndims = np.fromstring(next(sourceFile), dtype=int, sep=' ')
                    lineInd += 1
                else:
                    ndims = None
                lineInd += 2
                # TODO: Remove following line if Variant 2 is adopted
                # matIndToCheckEnd = len(matList) - 1
                for type, hLength in KNOWN_OCTAVE_TYPES.items():
                    if '# type: ' + type == typeLine[:-1]:
                        # TODO: Remove following lines when cleaning up the function
                        # for ind in range(hLength - 1):
                        #     firstNumericLine = next(sourceFile)
                        #     lineInd += 1
                        # if '(' in firstNumericLine and ')' in firstNumericLine:
                        #     numType = 'complex'
                        # else:
                        #     numType = 'real'
                        mat = {
                            'name': line[8:-1],
                            'type': type,
                            'startInd': startInd,
                            'headerLength': hLength,
                            'endInd': None,
                            'ndims': ndims,
                        }
                        matList.append(mat)
                # TODO: Remove following lines if Variant 2 is adopted
                # Variant 1
                # if matIndToCheckEnd >= 0 and matList[matIndToCheckEnd]['endInd'] is None:
                #     matList[matIndToCheckEnd]['endInd'] = startInd - 3
                #     matList[matIndToCheckEnd]['maxRows'] = (
                #        matList[matIndToCheckEnd]['endInd'] - matList[matIndToCheckEnd]['startInd']
                #         - matList[matIndToCheckEnd]['headerLength'] + 1
                #     )
            # Variant 2
            if line.strip() == '' and len(matList) > 0 and matList[-1]['endInd'] is None:
                matList[-1]['endInd'] = lineInd - 1
                matList[-1]['maxRows'] = matList[-1]['endInd'] - matList[-1]['startInd'] \
                    - matList[-1]['headerLength'] + 1
            # if line.strip() != '':
            #     lastLineWithContentInd = lineInd
        # # matList[-1]['endInd'] = totLines - 4 OLD
        # matList[-1]['endInd'] = lastLineWithContentInd
        # matList[-1]['maxRows'] = matList[-1]['endInd'] - matList[-1]['startInd'] \
        #     - matList[-1]['headerLength'] + 1
        totLines = lineInd
        for mat in matList:
            mat['skipRows'] = mat['startInd'] + mat['headerLength'] - 1
            mat['skipFooter'] = totLines - mat['endInd'] - 2
        # Correction to correctly read last variable in file with pd.read_csv(),
        # not clear why this is necesary
        matList[-1]['skipFooter'] = 0
    return matList


def load_octave_matrices(sourceFilePath, matNamesToLoad=None, colNames=None):
    matDefList = find_octave_all_matrices(sourceFilePath)
    matNamesAll = [mat['name'] for mat in matDefList]
    matStartInds = [mat['startInd'] for mat in matDefList]
    matHeaderLengths = [mat['headerLength'] for mat in matDefList]
    singleMatrixRequest = False
    if matNamesToLoad is None:
        matNamesToLoad = matNamesAll
    elif isinstance(matNamesToLoad, str):
        singleMatrixRequest = True
        matNamesToLoad = [matNamesToLoad]
    matList = {}
    for matDef in matDefList:
        if matDef['name'] in matNamesToLoad:
            if matDef['type'] in ['scalar', 'matrix']:
                matDf = pd.read_csv(
                    sourceFilePath, engine='python', delim_whitespace=True,
                    index_col=False, header=None, names=colNames,
                    skiprows=matDef['skipRows'], skipfooter=matDef['skipFooter']
                )
            elif matDef['type'] == 'complex matrix':
                matNp = np.loadtxt(
                    sourceFilePath, skiprows=matDef['skipRows'], max_rows=matDef['maxRows'],
                    dtype=np.complex128, converters={0: lambda s: np.fromstring(
                        s.decode("latin1").replace('(', '').replace(')', ''), sep=','
                    ).view(np.complex128)}
                )
                matDf = pd.DataFrame(data=matNp, columns=colNames)
            if colNames is not None:
                matList[matDef['name']] = matDf
            else:
                matNp = matDf.to_numpy().squeeze()
                if matDef['ndims'] is not None:
                    matNp = matNp.reshape(np.flip(matDef['ndims'])).transpose()
                elif matNp.shape == ():
                    matNp = matNp.item()
                matList[matDef['name']] = matNp
    if singleMatrixRequest:
        return matList[matNamesToLoad[0]]
    return matList
