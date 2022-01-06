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
# astraFilePath = './Elegant/FCCee_WP1p3/QuadOverRf/AstraReference/QuadOverRf.ini'
# astraFilePath = './Elegant/FCCee_WP1p3/QuadOverRf/AstraReference_VanishingLongEmit/QuadOverRf.ini'

zProjection = 500.   # [mm]
zCut = 500.   # [mm]
# astraFilePath = './Elegant/FCCee_WP1p3/QuadOverRf/AstraReference/QuadOverRf.0100.001'
# astraFilePath = './Elegant/FCCee_WP1p3/QuadOverRf/AstraReference_VanishingLongEmit/QuadOverRf.0100.001'
# astraFilePath = './Elegant/FCCee_WP1p3/QuadOverRf/AstraReference_NoRf/QuadOverRf.0100.001'
# astraFilePath = './Elegant/FCCee_WP1p3/QuadOverRf/AstraReference_NoRf_StrongerQuad4/QuadOverRf.0100.001'
astraFilePath = './Elegant/FCCee_WP1p3/QuadOverRf/AstraReference_StrongerQuadAndRf4/QuadOverRf.0100.001'
standardDf = bd.convert_astra_to_standard_df(
    astraFilePath, zProjection=zProjection, zCut=zCut,
    saveStandardFwf=True, verbose=True
)


#%%

# z0 = 0
# sddsFilepath = './Elegant/FCCee_WP1p3/QuadOverRf/Nslices1/QuadOverRf.bun'
# sddsFilepath = './Elegant/FCCee_WP1p3/QuadOverRf/Nslices1_VanishingLongEmit/QuadOverRf.bun'

# z0 = 2000.   # [mm]
# sddsFilepath = './Elegant/FCCee_WP1p3/QuadOverRf/PositiveDriftFirstOrder/QuadOverRf.out'
# sddsFilepath = './Elegant/FCCee_WP1p3/QuadOverRf/PositiveDriftSecondOrder/QuadOverRf.out'
# sddsFilepath = './Elegant/FCCee_WP1p3/QuadOverRf/PositiveDriftWithEmatrix/QuadOverRf.out'

z0 = 500.   # [mm]
# sddsFilepath = './Elegant/FCCee_WP1p3/QuadOverRf/NegativeDriftWithEmatrix/QuadOverRf.out'
# sddsFilepath = './Elegant/FCCee_WP1p3/QuadOverRf/Nslices1_NegativeDriftFirstOrder/QuadOverRf.out'
# sddsFilepath = './Elegant/FCCee_WP1p3/QuadOverRf/Nslices1_VanishingLongEmit/QuadOverRf.out'
# sddsFilepath = './Elegant/FCCee_WP1p3/QuadOverRf/Nslices1/QuadOverRf.out'
# sddsFilepath = './Elegant/FCCee_WP1p3/QuadOverRf/Nslices3/QuadOverRf.out'
# sddsFilepath = './Elegant/FCCee_WP1p3/QuadOverRf/Nslices6/QuadOverRf.out'
# sddsFilepath = './Elegant/FCCee_WP1p3/QuadOverRf/Nslices12/QuadOverRf.out'
# sddsFilepath = './Elegant/FCCee_WP1p3/QuadOverRf/Nslices100/QuadOverRf.out'
sddsFilepath = './Elegant/FCCee_WP1p3/QuadOverRf/Nslices1_NoRf/QuadOverRf_CrossDistr.out'
# sddsFilepath = './Elegant/FCCee_WP1p3/QuadOverRf/Nslices1_NoRf_Fringe/QuadOverRf_CrossDistr.out'
# sddsFilepath = './Elegant/FCCee_WP1p3/QuadOverRf/Nslices1_NoRf_StrongerQuad4/QuadOverRf.out'
# sddsFilepath = './Elegant/FCCee_WP1p3/QuadOverRf/Nslices1_StrongerQuadAndRf4/QuadOverRf.out'

standardDf = bd.convert_sdds_to_standard_df(
    sddsFilepath, z0=z0, pdgId=-11, Qbunch = 5.e-9, saveStandardFwf=True
)


#%%

fileBasePath = './Elegant/FCCee_WP1p3/QuadOverRf/AstraReference/Eabs_IdealTw'
freq = 2.9988e9
Lstructure = 0.5
zRes = 1.e-3
bd.generate_fieldmap_astra_ideal_tw(fileBasePath, freq, Lstructure, zRes)


#%%

elementName = 'QuadOverRf'
specLabel = ''
# specLabel = '_VanishingLongEmit'
# specLabel = '_NoRf'
# specLabel = '_StrongerQuadAndRf4'
L = 0.5   # [m]
Nslices = 1
enhancementFactor = 1.
quadGradient = enhancementFactor * 7.5   # [T/m]
quadOrder = 3
rfFreq = 2.9988e9   # [Hz]
rfPhase = -90.   # [deg]
rfVoltage = enhancementFactor * 9.   # [MV]
EkinIni = 200. # [MeV]
elegantInputFilePath = '/home/schaer_m/GIT_PSIPositronProduction/Elegant/FCCee_WP1p3/QuadOverRf/Nslices{0:d}{1:s}/QuadOverRf_Nslices{0:d}.lat'.format(
    Nslices, specLabel
)
bd.generate_lattice_quad_over_rf_elegant(
    elementName, L, Nslices, quadGradient, rfFreq, rfPhase, rfVoltage, EkinIni,
    elegantInputFilePath=elegantInputFilePath, quadOrder=quadOrder
)


#%%

xMax = 20.   # [mm]
yMax = 20.   # [mm]
p0 = 200.   # [MeV/c]
pzDelta = 100.   # [MeV/c]
outFilePath = './Elegant/FCCee_WP1p3/QuadOverRf/CrossDistribution'
standardDf = bd.generate_cross_distribution(
    xMax, yMax, p0, pzDelta, xPoints=7, yPoints=7, pzPoints=5,
    saveStandardFwf=True, outFilePath=outFilePath
)

astraDf = bd.convert_standard_df_to_astra(
    standardDf=standardDf, refParticleId=10,
    saveAstraDistr=True, outFilePath=outFilePath
)

bd.convert_standard_df_to_sdds(standardDf=standardDf, outFilePath=outFilePath)


#%%

# sdfFilePath = './Elegant/FCCee_WP1p3/QuadOverRf/AstraReference/QuadOverRf.ini.sdf_txt'
# sdfFilePath = './Elegant/FCCee_WP1p3/QuadOverRf/Nslices1/QuadOverRf.bun.sdf_txt'
sdfFilePath = '/afs/psi.ch/project/Pcubed/SimulationRuns/Geant4/000011/FCCeeTargetTracking_amor_leave_pdgId_-11.root.sdf_txt'

beam = bd.load_standard_fwf(sdfFilePath)
emitX = bd.compute_emittance(beam, 'x', correctOffsets=False)
print(emitX)


#%%

rootFilePath = '/afs/psi.ch/project/Pcubed/SimulationRuns/Geant4/000011/FCCeeTargetTracking.root'
bd.convert_fcceett_to_standard_df(
    rootFilePath, pdgId=-11, saveStandardFwf=True
)


#%%

octaveFilepath = './Elegant/FCCee_WP1p3/Distributions_Positrons_200MeV_Yongke/CTSB-N02-F100-E06-S0.5-T5.0_HTSTest_JNov04_SolC_CLICTW-Ztc200-Ri15-Bc0.50.dat'
z0 = 10e3   # [mm]
standardDf = bd.convert_octave_to_standard_df(
    octaveFilepath, z0=z0, pdgId=-11, Qbunch = 25.e-9, saveStandardFwf=True
)
