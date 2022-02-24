import RF_Track as rft
import BeamDynamics as bd
import numpy as np
import matplotlib.pyplot as plt


def create_solenoid(L, B0, Rin):   # m,T,m
    """Create solenoid fieldmap.

    Parameters
    ----------
    L : :obj:`float`
        Solenoid length in [m].
    Bz : :obj:`float`
        Solenoidal field in [T].
    Rin : :obj:`float`
        Solenoid internal radius (aperture) in [m].

    Returns
    -------
    solFieldmap : :obj:`RFTObj`
        RF-Track static magnetic fieldmap.

    """
    zAx = np.linspace(-2*L,2*L,200)
    lp = 0.5 * L + zAx   # [m]
    lm = 0.5 * L - zAx   # [m]
    Bax = 0.5*B0*(lm / np.hypot(Rin, lm) + lp / np.hypot(Rin, lp))   # [T]
    hz = zAx[1] - zAx[0]   # [m]
    solFieldmap = rft.Static_Magnetic_FieldMap_1d(Bax, hz)
    n_L = 95 / 0.236   # [turns/m]  (CLEAR solenoid)
    # T / mu0 / (1/m) = 795774.7150238460162655 A
    print(f'Solenoid current = {B0/n_L*795774.7150238460162655} A\n')
    return solFieldmap


# Bunch
Q = 4.   # [nC]
SIGMA_Z = 0.1   # [mm]
P_SPREAD = 1.   # [permil], momentum spread
P_REF = 200.   # [MeV/c]

Twiss = rft.Bunch6dT_twiss()
Twiss.emitt_x = 10e3   # [pi mm mrad]
Twiss.emitt_y = 10e3   # [pi mm mrad]
Twiss.beta_x = 1.   # [m]
Twiss.beta_y = 1.   # [m]
Twiss.alpha_x = 0.01
Twiss.alpha_y = 0.01
Twiss.sigma_z = SIGMA_Z   # [mm/c]
Twiss.sigma_d = P_SPREAD   # [permil]

N_PARTICLES = 10000   # Macro-particles
B0 = rft.Bunch6dT (
    rft.electronmass, Q*rft.nC, +1, P_REF, Twiss, N_PARTICLES
)
M0 = B0.get_phase_space()   # x [mm], Px [MeV/c], y [mm], Py [MeV/c], S [mm], Pz [MeV/c]
beamIn, _ = bd.convert_rftrack_to_standard_df(
    rftrackDf=M0, rftrackDfFormat='rftrack_Px_S', t=np.nan, pdgId=-11, Qbunch=Q*1e-9
)

# Solenoid
L_SOLENOID = 2.   # [m]
BZ_SOLENOID = 0.5   # [T]
R_IN_SOLENOID = 0.100   # [m]
solenoid = create_solenoid(L_SOLENOID, BZ_SOLENOID, R_IN_SOLENOID)
L_SPACING = L_SOLENOID + 0.3

# Drift
L_DRIFT = 10.    # [m]
EZ_ACC = 15e6   # [V/m]
drift = rft.Drift(L_DRIFT)
drift.set_static_Efield(0, 0, EZ_ACC)

# Tracking volume
vol = rft.Volume()
for solInd in range(0,3):
    vol.add(solenoid, 0, 0, solInd*L_SPACING, 'center')   # element, X, Y, Z in [m]
vol.add(drift, 0, 0, 0)

# Tracking options
TO = rft.TrackingOptions()
TO.dt_mm = 1.   # [mm/c]
TO.odeint_algorithm = 'rk2'   # 'rkf45', 'leapfrog', ...
TO.tt_dt_mm = 10.   # [mm/c], track the emittance every tt_dt_mm steps

# Plot Bz
fig1, ax1 = plt.subplots()
zAxis = np.linspace(vol.get_s0(), vol.get_s1(), 1000)   # [m]
Bz = []
for z in zAxis:
    [E, B] = vol.get_field(0, 0, z*1e3, 0)   # x,y,z,t (mm, mm/c)
    Bz.append(B[2])
Bz = np.array(Bz)
ax1.plot(zAxis, Bz)
ax1.set_xlabel('S [m]')
ax1.set_ylabel('Bz [T]')
ax1.grid()
plt.show()

# Tracking
B1 = vol.track(B0, TO)
M1 = B1.get_phase_space()   # x [mm], Px [MeV/c], y [mm], Py [MeV/c], S [mm], Pz [MeV/c]
# M1 = B1.get_phase_space('%X %xp %Y %yp %S %P')   # x [mm] Px [MeV/c]   y Py S Pz
TT = vol.get_transport_table('%mean_S %emitt_x %emitt_y %emitt_4d %sigma_X %sigma_Y')
beamOut, _ = bd.convert_rftrack_to_standard_df(
    rftrackDf=M1, rftrackDfFormat='rftrack_Px_S', t=np.nan, pdgId=-11, Qbunch=Q*1e-9
)

plotDefs = [
    {
        'varName1': 'x', 'varName2': 'y',
        'opacityHist': 0.6,
    },
    {
        'varName1': 'x', 'varName2': 'xp',
        'opacityHist': 0.6,
    },
    {
        'varName1': 'y', 'varName2': 'yp',
        'opacityHist': 0.6,
    },
    {
        'varName1': 'z', 'varName2': 'pz',
        'opacityHist': 0.6,
    },
    {
        'varName1': 'z', 'varName2': 'Ekin',
        'opacityHist': 0.6,
    },
]
ax = bd.plot_distr([beamIn, beamOut], plotDefs)
# # plt.ion()
plt.show()
# input("Press Enter to continue...")
# plt.close('all')

# TODO
# figure(4);
# clf ; hold on
# plot(TT(:,1)/1e3, TT(:,2));
# plot(TT(:,1)/1e3, TT(:,3));
# plot(TT(:,1)/1e3, TT(:,4));
# xlabel('S [m]');
# ylabel('emitt [mm.mrad]');
# legend('emitt x', 'emitt y', 'emitt 4d');
# figure(5);
# clf ; hold on
# plot(TT(:,1)/1e3, TT(:,5));
# plot(TT(:,1)/1e3, TT(:,6));
# xlabel('S [m]');
# ylabel('sigma [mm]');
# legend('sigma x', 'sigma y');
