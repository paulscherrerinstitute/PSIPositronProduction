#!/usr/bin/env python3


import sys
import BeamDynamics as bd


rootFileName = sys.argv[1]
#bd.convert_fcceett_to_standard_df(rootFileName, pdgId=-11, saveStandardFwf=True)
bd.convert_fcceett_to_standard_df(rootFileName, saveStandardFwf=True)
