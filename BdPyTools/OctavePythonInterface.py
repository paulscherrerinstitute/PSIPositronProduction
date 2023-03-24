import numpy as np
import pandas as pd
import timeit


def find_octave_all_matrices(sourceFilePath):
    KNOWN_OCTAVE_TYPES = {
        # typeName: headerLength
        'scalar': 2,
        'matrix': 4,
        'complex matrix': 4,
    }
    with open(sourceFilePath, 'r') as sourceFile:
        matList = []
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
                for type, hLength in KNOWN_OCTAVE_TYPES.items():
                    if '# type: ' + type == typeLine[:-1]:
                        mat = {
                            'name': line[8:-1],
                            'type': type,
                            'startInd': startInd,
                            'headerLength': hLength,
                            'endInd': None,
                            'ndims': ndims,
                        }
                        matList.append(mat)
            if line.strip() == '' and len(matList) > 0 and matList[-1]['endInd'] is None:
                matList[-1]['endInd'] = lineInd - 1
                matList[-1]['maxRows'] = matList[-1]['endInd'] - matList[-1]['startInd'] \
                    - matList[-1]['headerLength'] + 1
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
                converters = None
                dtype = np.float64
            elif matDef['type'] == 'complex matrix':
                dtype = np.complex128

                def string_to_complex(byteStr):
                    charStr = byteStr.decode('latin1')[1:-1]
                    complexArray = np.fromstring(charStr, sep=',').view(np.complex128)
                    return complexArray
                converters = string_to_complex
            # execTimeStart = timeit.default_timer()
            matNp = np.loadtxt(
                sourceFilePath, skiprows=matDef['skipRows'], max_rows=matDef['maxRows'],
                dtype=dtype, converters=converters)
            # execTimeStop = timeit.default_timer()
            # print('Time for np.loadtxt(): ', execTimeStop-execTimeStart, 's')
            if colNames is not None:
                matList[matDef['name']] = pd.DataFrame(data=matNp, columns=colNames)
            else:
                if matDef['ndims'] is not None:
                    matNp = matNp.reshape(np.flip(matDef['ndims'])).transpose()
                elif matNp.shape == ():
                    matNp = matNp.item()
                matList[matDef['name']] = matNp
    if singleMatrixRequest:
        return matList[matNamesToLoad[0]]
    return matList
