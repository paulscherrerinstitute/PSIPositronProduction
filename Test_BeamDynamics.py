import BeamDynamics as bd


# bd.convert_irina_distr_to_standard_df(
#     'Geant4/FcceeTarget_StartingExample/run_Injector/ex_gen1.dat',
#     saveStandardCsv=True
# )

# sddsFilepath = './Elegant/FCCee_WP1p3/QuadOverRf1/QuadOverRf1.bun'
# sddsFilepath = './Elegant/FCCee_WP1p3/QuadOverRf1/QuadOverRf1.out'
sddsFilepath = './Elegant/FCCee_WP1p3/QuadOverRf1/PositiveDrift/QuadOverRf1.out'
# sddsFilepath = './Elegant/FCCee_WP1p3/QuadOverRf1/PositiveDriftWithEmatrix/QuadOverRf1.out'
# sddsFilepath = './Elegant/FCCee_WP1p3/QuadOverRf1/NegativeDriftWithEmatrix/QuadOverRf1.out'
standardDf = bd.convert_sdds_to_standard_df(
    sddsFilepath, pdgId = -11, z = 0, Qbunch = 5.e-9, saveStandardFwf=True
)
