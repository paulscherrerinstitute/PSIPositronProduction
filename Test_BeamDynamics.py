import BeamDynamics as bd
from importlib import reload
reload(bd)


#%%

# bd.convert_irina_distr_to_standard_df(
#     'Geant4/FcceeTarget_StartingExample/run_Injector/ex_gen1.dat',
#     saveStandardCsv=True
# )


#%%

# zProjection = None
# zCut = None
# astraFilePath = './Elegant/FCCee_WP1p3/QuadOverRf1/AstraReference/QuadOverRf.ini'

zProjection = 500.   # [mm]
zCut = 500.   # [mm]
# astraFilePath = './Elegant/FCCee_WP1p3/QuadOverRf1/AstraReference/QuadOverRf.0100.001'
# astraFilePath = './Elegant/FCCee_WP1p3/QuadOverRf1/AstraReference_VanishingLongEmit/QuadOverRf.ini'
# astraFilePath = './Elegant/FCCee_WP1p3/QuadOverRf1/AstraReference_VanishingLongEmit/QuadOverRf1.0050.001'
# astraFilePath = './Elegant/FCCee_WP1p3/QuadOverRf1/AstraReference_LongEmit/QuadOverRf.ini'
# astraFilePath = './Elegant/FCCee_WP1p3/QuadOverRf1/AstraReference_LongEmit/QuadOverRf.0050.001'
astraFilePath = './Elegant/FCCee_WP1p3/QuadOverRf1/AstraReference_NoRf/QuadOverRf.0100.001'
# astraFilePath = './Elegant/FCCee_WP1p3/QuadOverRf1/AstraReference_StrongerQuad4/QuadOverRf.0050.001'
# astraFilePath = './Elegant/FCCee_WP1p3/QuadOverRf1/AstraReference_StrongerQuadAndRf4/QuadOverRf.0050.001'
# astraFilePath = './Elegant/FCCee_WP1p3/QuadOverRf1/AstraReference_OnlyQuad/QuadOverRf.0100.001'
bd.convert_astra_to_standard_df(
    astraFilePath, zProjection=zProjection, zCut=zCut,
    saveStandardFwf=True, verbose=True
)


#%%

# z0 = 0
# sddsFilepath = './Elegant/FCCee_WP1p3/QuadOverRf1/Nslices1/QuadOverRf.bun'
# sddsFilepath = './Elegant/FCCee_WP1p3/QuadOverRf1/Nslices1_LongEmit/QuadOverRf.bun'
# sddsFilepath = './Elegant/FCCee_WP1p3/QuadOverRf1/Nslices1_VanishingLongEmit/QuadOverRf.bun'

z0 = 500.   # [mm]
# sddsFilepath = './Elegant/FCCee_WP1p3/QuadOverRf1/PositiveDrift/QuadOverRf.out'
# sddsFilepath = './Elegant/FCCee_WP1p3/QuadOverRf1/PositiveDriftWithEmatrix/QuadOverRf.out'
# sddsFilepath = './Elegant/FCCee_WP1p3/QuadOverRf1/NegativeDriftWithEmatrix/QuadOverRf.out'
# sddsFilepath = './Elegant/FCCee_WP1p3/QuadOverRf1/Nslices1_VanishingLongEmit/QuadOverRf.out'
# sddsFilepath = './Elegant/FCCee_WP1p3/QuadOverRf1/Nslices1/QuadOverRf.out'
# sddsFilepath = './Elegant/FCCee_WP1p3/QuadOverRf1/Nslices1_LongEmit/QuadOverRf.out'
# sddsFilepath = './Elegant/FCCee_WP1p3/QuadOverRf1/Nslices1/QuadOverRf.out'
# sddsFilepath = './Elegant/FCCee_WP1p3/QuadOverRf1/Nslices3/QuadOverRf.out'
# sddsFilepath = './Elegant/FCCee_WP1p3/QuadOverRf1/Nslices6/QuadOverRf.out'
# sddsFilepath = './Elegant/FCCee_WP1p3/QuadOverRf1/Nslices12/QuadOverRf.out'
# sddsFilepath = './Elegant/FCCee_WP1p3/QuadOverRf1/Nslices100/QuadOverRf.out'
sddsFilepath = './Elegant/FCCee_WP1p3/QuadOverRf1/Nslices1_NoRf/QuadOverRf.out'
# sddsFilepath = './Elegant/FCCee_WP1p3/QuadOverRf1/Nslices1_NoRf_StrongerQuad4/QuadOverRf.out'
# sddsFilepath = './Elegant/FCCee_WP1p3/QuadOverRf1/Nslices1_StrongerQuad4/QuadOverRf.out'
# sddsFilepath = './Elegant/FCCee_WP1p3/QuadOverRf1/Nslices1_StrongerQuadAndRf4/QuadOverRf.out'
# sddsFilepath = './Elegant/FCCee_WP1p3/QuadOverRf1/Nslices100_StrongerQuadAndRf4/QuadOverRf.out'

standardDf = bd.convert_sdds_to_standard_df(
    sddsFilepath, z0=z0, pdgId=-11, Qbunch = 5.e-9, saveStandardFwf=True
)


#%%

fileBasePath = './Elegant/FCCee_WP1p3/QuadOverRf1/AstraReference/Eabs_IdealTw'
freq = 2.9988e9
Lstructure = 0.5
zRes = 1.e-3
bd.generate_fieldmap_astra_ideal_tw(fileBasePath, freq, Lstructure, zRes)


#%%

elementName = 'QuadOverRf'
# specLabel = ''
# specLabel = '_VanishingLongEmit'
specLabel = '_NoRf'
L = 0.5   # [m]
Nslices = 1
quadGradient = 1. * 7.5   # [T/m]
rfFreq = 2.9988e9   # [Hz]
rfPhase = -90.   # [deg]
rfVoltage = 1. * 9.   # [MV]
EkinIni = 200. # [MeV]
elegantInputFilePath = '/home/schaer_m/GIT_PSIPositronProduction/Elegant/FCCee_WP1p3/QuadOverRf1/Nslices{0:d}{1:s}/QuadOverRf_Nslices{0:d}.lat'.format(
    Nslices, specLabel
)
bd.generate_lattice_quad_over_rf_elegant(
    elementName, L, Nslices, quadGradient, rfFreq, rfPhase, rfVoltage, EkinIni,
    elegantInputFilePath=elegantInputFilePath
)


#%%

xMax = 10.   # [mm]
yMax = 10.   # [mm]
p0 = 200.   # [MeV/c]
pzDelta = 100.   # [MeV/c]
outFilePath = './Elegant/FCCee_WP1p3/QuadOverRf1/CrossDistribution'
standardDf = bd.generate_cross_distribution(
    xMax, yMax, p0, pzDelta, xPoints=5, yPoints=5, pzPoints=5,
    saveStandardFwf=True, outFilePath=outFilePath
)

bd.convert_standard_df_to_astra(
    standardDf=standardDf, refParticleId=10,
    saveAstraDist=True, outFilePath=outFilePath
)

