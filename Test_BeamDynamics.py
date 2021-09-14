import BeamDynamics as bd


# bd.convert_irina_distr_to_standard_df(
#     'Geant4/FcceeTarget_StartingExample/run_Injector/ex_gen1.dat',
#     saveStandardCsv=True
# )

standardDf = bd.convert_sdds_to_standard_df(
    './Elegant/FCCee_WP1p3/QuadOverRf_BasicExample_1/QuadOverRf1.out',
    pdgId = -11, z = 0, Qbunch = 5.e-9, saveStandardFwf=True
)
