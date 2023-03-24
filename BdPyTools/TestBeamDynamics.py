"""
Unit Testing of Module Beam Dynamics
====================================

"""

import unittest
import os
import pathlib as pl
import OctavePythonInterface as opi
import BeamDynamics as bd


def assert_file_generated(outFilePath):
    """Assert if a file has been generated or not.
    Typically used to test format conversions.

    Parameters
    ----------
    outFilePath : :obj:`str`
        Path of the file to be searched for.

    """
    outFile = pl.Path(outFilePath)
    # TODO: Add resolve below? outFile.resolve().is_file()
    if not outFile.is_file():
        raise AssertionError(f'File {outFilePath} does not exist.')
    os.remove(outFilePath)


class TestBeamDynamics(unittest.TestCase):
    """Unit testing for the class BeamDynamics."""

    # @classmethod
    # def setUpClass(cls):
    #     cls.testData = loadourDataFunc(...)

    def test_astra_to_standard_1(self):
        """Test conversion from ASTRA format to standard data frame, not discarding any particle."""
        discardLostParticles = True
        zProjection = None
        zCut = None
        astraFilePath = 'Data/Astra/RUN_1911_164058.1380.001'
        standardDf, outFwfPath = bd.convert_astra_to_standard_df(
            astraFilePath, discardLostParticles=discardLostParticles, zProjection=zProjection,
            zCut=zCut, outFwfPath=pl.Path(astraFilePath).stem, verbose=True
        )
        assert_file_generated(outFwfPath)
        self.assertEqual(standardDf.shape, (14558, 14))
        self.assertAlmostEqual(standardDf.loc[0, 'px'], 0.0)
        self.assertAlmostEqual(standardDf.loc[30911, 'py'], -1.726)
        self.assertAlmostEqual(standardDf.loc[7954, 'pz'], 57.33000000000001)

    def test_standard_to_astra_1(self):
        """Test conversion from standard data frame to ASTRA format."""
        sourceFilePath = 'Data/RFTrack/CaptureLinac/CaptureLinacUpTo200MeV_LBandLargeR' + \
            '_RealSolenoids_Type1and2_TargetAt35mm/DistrOut_After1stTracking_6d.sdf_txt'
        astraFilePath = 'DistrOut_After1stTracking_6d.0001.001'
        df = bd.load_standard_fwf(sourceFilePath)
        totPositronCharge = 5e-9  # [C]
        df['Q'] = totPositronCharge / df.shape[0]
        bd.convert_standard_df_to_astra(standardDf=df, outFilePath=pl.Path(astraFilePath).stem[:-5])
        assert_file_generated(astraFilePath)

    def test_rftrack_to_standard_df_1(self):
        """Test conversion from RF-Track format to standard data frame, example 1."""
        sourceFilePath = 'Data/RFTrack/PositronsAt200MeV_Yongke/' + \
            'CTSB-N02-F100-E06-S0.5-T5.0_HTSTest_JNov04_SolC_PSISW-Ztc200-Ri15-Bc1.50.dat'
        octaveMatrixName = 'A_RF'
        # TODO: Following 2 values are very approximative
        s = 10e3   # [mm]
        Qbunch = 25.e-9   # [C]
        standardDf, outFwfPath = bd.convert_rftrack_to_standard_df(
            sourceFilePath=sourceFilePath, sourceFormat='octave', octaveMatrixName=octaveMatrixName,
            rftrackDfFormat='rftrack_xp_t', s=s, pdgId=-11, Qbunch=Qbunch,
            outFwfPath=pl.Path(sourceFilePath).stem
        )
        assert_file_generated(outFwfPath)
        self.assertAlmostEqual(standardDf.loc[0, 'px'], 2.286897262e-01)
        self.assertAlmostEqual(standardDf.loc[54324, 'py'], -1.931028202e+00)

    def test_standard_to_rftrack(self):
        """Test conversion from standard data frame to Octave RF-Track format."""
        sdfFilePath = 'Data/RFTrack/CaptureLinac/' + \
            'CaptureLinacUpTo200MeV_LBandLargeR_RealSolenoids_Type1and2/' + \
            'DistrOut_After1stTracking_6d.sdf_txt'
        rftFilePath = 'DistrOut_After1stTracking_6d.dat'
        rftrackDfFormat = 'rftrack_xp_t'
        bd.convert_standard_df_to_rftrack(
            sourceFilePath=sdfFilePath, rftrackDfFormat=rftrackDfFormat,
            outFilePath=pl.Path(rftFilePath).stem
        )
        assert_file_generated(rftFilePath)

    def test_rftrack_back_and_forth(self):
        """Test conversion from RF-Track format to standard data frame and back."""
        # 'CTSB-N02-F100-E06-S0.5-T5.0_HTSTest_JNov04_SolC_CLICTW-Ztc200-Ri15-Bc0.50.dat'
        # 'CTSB-N02-F100-E06-S0.5-T5.0_FCTest_Pavel_SolC_CLICTW-OptionA08B7-Bc0.50.dat'
        octaveMatrixName = 'A_RF'
        rftrackDfFormat = 'rftrack_xp_t'
        colNames = bd.FILE_TYPES_SPECS[rftrackDfFormat]['columnOrder']
        sourceFilePath = 'Data/RFTrack/PositronsAt200MeV_Yongke/' + \
            'CTSB-N02-F100-E06-S0.5-T5.0_HTSTest_JNov04_SolC_PSISW-Ztc200-Ri15-Bc1.50.dat'
        indsAnglesNotConformal = [88346, 93448]
        # indPzNegative = [7666, 22069, 88346, 93448]
        rftDfOriginal = opi.load_octave_matrices(
            sourceFilePath, matNamesToLoad='A_RF', colNames=colNames
        )
        # TODO: Following 2 values are very approximative
        s = 10e3   # [mm]
        Qbunch = 25.e-9   # [C]
        standardDf, _ = bd.convert_rftrack_to_standard_df(
            sourceFilePath=sourceFilePath, octaveMatrixName=octaveMatrixName,
            rftrackDfFormat=rftrackDfFormat, s=s, pdgId=-11, Qbunch=Qbunch, outFwfPath=None
        )
        # TODO check t= vs z= in rftrack script
        # with open(sd.build_data_path(REL_PATH, 'filterSpces.json'), 'r') as filterSpecsFile:
        #     filterSpecsList = json.load(filterSpecsFile)
        # filterSpecs = filterSpecsList[FILTER_SPECS_SELECTOR]['filterSpecs']
        # beamIn = bd.filter_distr(beamIn, filterSpecs)
        rftDf, _ = bd.convert_standard_df_to_rftrack(
            standardDf=standardDf, rftrackDfFormat=rftrackDfFormat
        )
        numDecimalPlaces = 6
        dfDiff = rftDf.drop(indsAnglesNotConformal).round(numDecimalPlaces) \
            .compare(rftDfOriginal.drop(indsAnglesNotConformal).round(numDecimalPlaces))
        self.assertEqual(dfDiff.shape[0], 0)

    def test_standard_df_to_sdds(self):
        """Test conversion from standard data frame to SDDS format."""
        sdfFilePath = 'Data/RFTrack/CaptureLinac/' + \
            'PositronLinac_15RFStruct_LBandLargeR_SolenoidsType1and2/' + \
            'DistrOut_After2ndTracking_6d.sdf_txt'
        sddsFilePath = 'DistrOut_After2ndTracking_6d.sdds'
        bd.convert_standard_df_to_sdds(
            sourceFilePath=sdfFilePath, outFilePath=pl.Path(sddsFilePath).stem
        )
        assert_file_generated(sddsFilePath)

    def test_plot_1(self):
        """Test plotting of phase space, standard example."""
        discardLostParticles = True
        zProjection = 1380.
        zCut = None
        astraFilePath = 'Data/Astra/RUN_1911_164058.1380.001'
        standardDf, _ = bd.convert_astra_to_standard_df(
            astraFilePath, discardLostParticles=discardLostParticles,
            zProjection=zProjection, zCut=zCut
        )
        plotDefs = bd.set_plot_defs_from_distrs(
            [standardDf], setNames=['TransvPlane', 'TransvPsAngles', 'LongPsT']
        )
        ax = bd.plot_distr([standardDf], plotDefs)

    def test_plot_2(self):
        """Test plotting of phase space,
        example with point-like distributions on certain coordinates."""
        sourceFilePath = 'Data/RFTrack/CaptureLinac/' + \
            'CaptureLinacUpTo200MeV_LBandLargeR_RealSolenoids_Type1and2_TargetAt35mm/' + \
            'DistrOut_PositronLinac.dat'
        standardDf, _ = bd.convert_rftrack_to_standard_df(
            sourceFilePath=sourceFilePath, sourceFormat='octave', octaveMatrixName='A_PL',
            rftrackDfFormat='rftrack_xp_t', pdgId=-11
        )
        plotDefs = bd.set_plot_defs_from_distrs(
            [standardDf], setNames=['TransvPlane', 'TransvPsAngles', 'LongPsT']
        )
        ax = bd.plot_distr([standardDf], plotDefs)

    # zProjection = 500.   # [mm]
    # zCut = 500.   # [mm]
    # # astraFilePath = './Elegant/FCCee_WP1p3/QuadOverRf/AstraReference/QuadOverRf.0100.001'
    # # astraFilePath = './Elegant/FCCee_WP1p3/QuadOverRf/AstraReference_VanishingLongEmit/QuadOverRf.0100.001'
    # # astraFilePath = './Elegant/FCCee_WP1p3/QuadOverRf/AstraReference_NoRf/QuadOverRf.0100.001'
    # # astraFilePath = './Elegant/FCCee_WP1p3/QuadOverRf/AstraReference_NoRf_StrongerQuad4/QuadOverRf.0100.001'
    # astraFilePath = './Elegant/FCCee_WP1p3/QuadOverRf/AstraReference_StrongerQuadAndRf4/QuadOverRf.0100.001'
    # standardDf = bd.convert_astra_to_standard_df(
    #     astraFilePath, zProjection=zProjection, zCut=zCut,
    #     saveStandardFwf=True, verbose=True
    # )

    #%%
    #
    # bd.convert_irina_distr_to_standard_df(
    #     'Geant4/FcceeTarget_StartingExample/run_Injector/ex_gen1.dat',
    #     saveStandardCsv=True
    # )

    #%%

    # # z0 = 0
    # # sddsFilepath = './Elegant/FCCee_WP1p3/QuadOverRf/Nslices1/QuadOverRf.bun'
    # # sddsFilepath = './Elegant/FCCee_WP1p3/QuadOverRf/Nslices1_VanishingLongEmit/QuadOverRf.bun'

    # # z0 = 2000.   # [mm]
    # # sddsFilepath = './Elegant/FCCee_WP1p3/QuadOverRf/PositiveDriftFirstOrder/QuadOverRf.out'
    # # sddsFilepath = './Elegant/FCCee_WP1p3/QuadOverRf/PositiveDriftSecondOrder/QuadOverRf.out'
    # # sddsFilepath = './Elegant/FCCee_WP1p3/QuadOverRf/PositiveDriftWithEmatrix/QuadOverRf.out'

    # z0 = 500.   # [mm]
    # # sddsFilepath = './Elegant/FCCee_WP1p3/QuadOverRf/NegativeDriftWithEmatrix/QuadOverRf.out'
    # # sddsFilepath = './Elegant/FCCee_WP1p3/QuadOverRf/Nslices1_NegativeDriftFirstOrder/QuadOverRf.out'
    # # sddsFilepath = './Elegant/FCCee_WP1p3/QuadOverRf/Nslices1_VanishingLongEmit/QuadOverRf.out'
    # # sddsFilepath = './Elegant/FCCee_WP1p3/QuadOverRf/Nslices1/QuadOverRf.out'
    # # sddsFilepath = './Elegant/FCCee_WP1p3/QuadOverRf/Nslices3/QuadOverRf.out'
    # # sddsFilepath = './Elegant/FCCee_WP1p3/QuadOverRf/Nslices6/QuadOverRf.out'
    # # sddsFilepath = './Elegant/FCCee_WP1p3/QuadOverRf/Nslices12/QuadOverRf.out'
    # # sddsFilepath = './Elegant/FCCee_WP1p3/QuadOverRf/Nslices100/QuadOverRf.out'
    # sddsFilepath = './Elegant/FCCee_WP1p3/QuadOverRf/Nslices1_NoRf/QuadOverRf_CrossDistr.out'
    # # sddsFilepath = './Elegant/FCCee_WP1p3/QuadOverRf/Nslices1_NoRf_Fringe/QuadOverRf_CrossDistr.out'
    # # sddsFilepath = './Elegant/FCCee_WP1p3/QuadOverRf/Nslices1_NoRf_StrongerQuad4/QuadOverRf.out'
    # # sddsFilepath = './Elegant/FCCee_WP1p3/QuadOverRf/Nslices1_StrongerQuadAndRf4/QuadOverRf.out'

    # standardDf = bd.convert_sdds_to_standard_df(
    #     sddsFilepath, z0=z0, pdgId=-11, Qbunch = 5.e-9, saveStandardFwf=True
    # )

    #%%
    #
    # fileBasePath = './Elegant/FCCee_WP1p3/QuadOverRf/AstraReference/Eabs_IdealTw'
    # freq = 2.9988e9
    # Lstructure = 0.5
    # zRes = 1.e-3
    # bd.generate_fieldmap_astra_ideal_tw(fileBasePath, freq, Lstructure, zRes)

    #%%
    #
    # elementName = 'QuadOverRf'
    # specLabel = ''
    # # specLabel = '_VanishingLongEmit'
    # # specLabel = '_NoRf'
    # # specLabel = '_StrongerQuadAndRf4'
    # L = 0.5   # [m]
    # Nslices = 1
    # enhancementFactor = 1.
    # quadGradient = enhancementFactor * 7.5   # [T/m]
    # quadOrder = 3
    # rfFreq = 2.9988e9   # [Hz]
    # rfPhase = -90.   # [deg]
    # rfVoltage = enhancementFactor * 9.   # [MV]
    # EkinIni = 200. # [MeV]
    # elegantInputFilePath = '/home/schaer_m/GIT_PSIPositronProduction/Elegant/FCCee_WP1p3/QuadOverRf/Nslices{0:d}{1:s}/QuadOverRf_Nslices{0:d}.lat'.format(
    #     Nslices, specLabel
    # )
    # bd.generate_lattice_quad_over_rf_elegant(
    #     elementName, L, Nslices, quadGradient, rfFreq, rfPhase, rfVoltage, EkinIni,
    #     elegantInputFilePath=elegantInputFilePath, quadOrder=quadOrder
    # )

    #%%
    #
    # xMax = 20.   # [mm]
    # yMax = 20.   # [mm]
    # p0 = 200.   # [MeV/c]
    # pzDelta = 100.   # [MeV/c]
    # outFilePath = './Elegant/FCCee_WP1p3/QuadOverRf/CrossDistribution'
    # standardDf = bd.generate_cross_distribution(
    #     xMax, yMax, p0, pzDelta, xPoints=7, yPoints=7, pzPoints=5,
    #     saveStandardFwf=True, outFilePath=outFilePath
    # )

    # bd.convert_standard_df_to_sdds(standardDf=standardDf, outFilePath=outFilePath)

    #%%
    #
    # # sdfFilePath = './Elegant/FCCee_WP1p3/QuadOverRf/AstraReference/QuadOverRf.ini.sdf_txt'
    # # sdfFilePath = './Elegant/FCCee_WP1p3/QuadOverRf/Nslices1/QuadOverRf.bun.sdf_txt'
    # sdfFilePath = '/afs/psi.ch/project/Pcubed/SimulationRuns/Geant4/000011/FCCeeTargetTracking_amor_leave_pdgId_-11.root.sdf_txt'

    # beam = bd.load_standard_fwf(sdfFilePath)
    # emitX = bd.compute_emittance(beam, 'x', correctOffsets=False)
    # print(emitX)

    #%%
    #
    # rootFilePath = '/afs/psi.ch/project/Pcubed/SimulationRuns/Geant4/000011/FCCeeTargetTracking.root'
    # bd.convert_fcceett_to_standard_df(
    #     rootFilePath, pdgId=-11, saveStandardFwf=True
    # )

    #%%
    #
    # sourceFilePath = '../FCCeeInjectorBeamApp/BeamDistrs/Positrons_200MeV_Yongke/CTSB-N02-F100-E06-S0.5-T5.0_HTSTest_JNov04_SolC_CLICTW-Ztc200-Ri15-Bc0.50.dat.sdf_txt'
    # filterSpecs = {
    #     'x': (-25., 25.),
    #     'xp': (-30., 30.),
    #     'y': (-25., 25.),
    #     'yp': (-30., 30.),
    #     't': (62.8, 63.1),
    #     'pz': (50., 400.),
    # }
    # astraRefParticle = {
    #     'x': 0.,
    #     'px': 0.,
    #     'y': 0.,
    #     'py': 0.,
    #     'pz': 200.,
    #     't': 62.87,
    #     'pdgId': -11,
    # }
    #
    # # sourceFilePath = '../FCCeeInjectorBeamApp/BeamDistrs/Positrons_200MeV_Yongke/CTSB-N02-F100-E06-S0.5-T5.0_HTSTest_JNov04_SolC_PSISW-Ztc200-Ri15-Bc1.50.dat.sdf_txt'
    # # filterSpecs = {
    # #     'x': (-25., 25.),
    # #     'xp': (-30., 30.),
    # #     'y': (-25., 25.),
    # #     'yp': (-30., 30.),
    # #     't': (58.6, 58.9),
    # #     'pz': (0, 400.),
    # # }
    # # astraRefParticle = {
    # #     'x': 0.,
    # #     'px': 0.,
    # #     'y': 0.,
    # #     'py': 0.,
    # #     'pz': 225.,
    # #     't': 58.73,
    # #     'pdgId': -11,
    # # }
    #
    # # bd.convert_standard_df_to_placet(sourceFilePath=sourceFilePath, savePlacetDistr=True)
    #
    # standardDf = bd.load_standard_fwf(sourceFilePath)
    # standardDf = bd.filter_distr(standardDf, filterSpecs)
    #
    # # bd.convert_standard_df_to_placet(standardDf=standardDf, savePlacetDistr=True)
    #
    # astraRefParticle['z'] = standardDf['z'][0]
    # astraRefParticle['Q'] = standardDf['Q'][0]
    # standardDf = standardDf.append(astraRefParticle, ignore_index=True)
    # # bd.convert_standard_df_to_astra(
    # #     standardDf=standardDf, refParticleId=standardDf.shape[0]-1,
    # #     saveAstraDistr=True, outFilePath=sourceFilePath
    # # )
    #
    # bd.convert_standard_df_to_sdds(
    #     standardDf=standardDf, refParticleId=standardDf.shape[0]-1,
    #     outFilePath=sourceFilePath
    # )

    #%%
    # # Astra to Placet
    #
    # sourceFilePath = '/home/tia/tmp/EmitGrowthInDriftSpace/NicoDistr_Z20p57m/SingleBucket/RUN_2501_121416.2057_SingleBucket.001'
    # zProjection = 20.57e3   # [mm]
    #
    # standardDf = bd.convert_astra_to_standard_df(
    #     sourceFilePath, zProjection=zProjection,
    #     saveStandardFwf=False, verbose=True
    # )
    # bd.convert_standard_df_to_placet(
    #     standardDf=standardDf, outFilePath=sourceFilePath
    # )


if __name__ == '__main__':
    unittest.main()
