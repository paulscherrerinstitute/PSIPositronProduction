"""
Simulation Data
====================================

Definitions and tools to retrieve simulation data.

"""


import os


DATA_BASE_PATH = '/afs/psi.ch/project/Pcubed'


def build_data_path(*relPathPieces):
    """Build absolute data path
    starting from input relative paths and module constants.

    Parameters
    ----------
    *relPathPieces : :obj:`str`
        One or more relative paths that will be concatenated.
        Cannot start with a separator.

    Returns
    -------
    absPath : :obj:`str`
        Absolute path.

    """
    for relPathPiece in relPathPieces:
        if relPathPiece[0] in ['/', '\\']:
            raise ValueError(
                f'relPathPiece cannot start with a path separator ({relPathPiece[0]} found).'
            )
    absPath = os.path.join(DATA_BASE_PATH, *relPathPieces)
    return absPath
