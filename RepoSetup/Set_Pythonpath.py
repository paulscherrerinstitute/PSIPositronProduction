#!/usr/bin/env python3



import os


pythonEnvVariables = (
    'PYTHONPATH',
    'JUPYTER_PATH',
)

moduleRelPaths = (
    '.',
    './trial',
)

repoPath = os.getcwd()

for pev in pythonEnvVariables:
    if pev not in os.environ:
        firstEntry = True
    else:
        firstEntry = False
    for mrp in moduleRelPaths:
        moduleAbsPath = os.path.join(repoPath, mrp)
        if firstEntry:
            os.environ[pev] = moduleAbsPath
            firstEntry = False
        elif mrp not in os.environ[pev]:
            os.environ[pev] += os.pathsep + moduleAbsPath
    os.system('export '+pev+'='+os.environ[pev])
    print('Environment variable '+pev+' is:\n'+os.popen('echo $'+pev).read())
