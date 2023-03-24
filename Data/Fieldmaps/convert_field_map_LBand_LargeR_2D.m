clear all;

% filenameEcomplex = 'pLinacF3_full44cells_Ecomp_YZplane_dy2mm_dz0p1L.fld';
% filenameHcomplex = 'pLinacF3_full44cells_Hcomp_YZplane_dy2mm_dz0p1L.fld';
% matShape = [461, 16];
% filenameOctaveFieldMap = 'pLinacF3_full44cells_YZplane_dy2mm_dz0p1L.dat';
% or
filenameEcomplex = 'pLinacF3_full44cells_m2_Ecomp_YZplane_dy1mm_dz0p025L.fld';
filenameHcomplex = 'pLinacF3_full44cells_m2_Hcomp_YZplane_dy1mm_dz0p025L.fld';
matShape = [1841, 31];
filenameOctaveFieldMap = 'pLinacF3_full44cells_m2_YZplane_dy1mm_dz0p025L.dat';

fid = fopen (filenameEcomplex, 'r');
rawE = textscan(fid, '%f%f%f%f%f%f%f%f%f', 'HeaderLines', 2);
fclose(fid);
rawE = cell2mat(rawE);
size(rawE)
Ex = rawE(:, 4) + j * rawE(:, 5);
Ey = rawE(:, 6) + j * rawE(:, 7);
Ez = rawE(:, 8) + j * rawE(:, 9);

fid = fopen (filenameHcomplex, 'r');
rawH = textscan(fid, '%f%f%f%f%f%f%f%f%f', 'HeaderLines', 2);
fclose(fid);
rawH = cell2mat(rawH);
size(rawH)
Hx = rawH(:, 4) + j * rawH(:, 5);
Hy = rawH(:, 6) + j * rawH(:, 7);
Hz = rawH(:, 8) + j * rawH(:, 9);

for axInd = 1:3
    if sum(rawE(:, axInd) ~= rawH(:, axInd)) > 0
        error('Meshes do not correspond!')
    endif
endfor
X = rawE(:, 1);
Y = rawE(:, 2);
Z = rawE(:, 3);

amplitudeScaling = 3.3513e7 / 5.1797e3;
globalPhaseAlignment = exp(pi*j);
mu0 = 1.25663706212e-6;  % [H/m]
THETA = reshape(X, matShape);
R = reshape(Y, matShape);
Z = reshape(Z, matShape) + 0.11942;
Etheta = reshape(Ex, matShape) * amplitudeScaling * globalPhaseAlignment;
Er = reshape(Ey, matShape) * amplitudeScaling * globalPhaseAlignment;
Ez = reshape(Ez, matShape) * amplitudeScaling * globalPhaseAlignment;
Btheta = reshape(Hx, matShape) * amplitudeScaling * globalPhaseAlignment * mu0;
Br = reshape(Hy, matShape) * amplitudeScaling * globalPhaseAlignment * mu0;
Bz = reshape(Hz, matShape) * amplitudeScaling * globalPhaseAlignment * mu0;
size(R)
size(THETA)
size(Z)
size(Er)
size(Etheta)
size(Ez)
size(Br)
size(Btheta)
size(Bz)

frequency = 2e9  % [Hz]
phase_advance = 9./10. * pi  % [rad]
wave_direction = +1
cell_length = 2.99792458e8 / frequency * 9./20.  % [m]
tot_cells = 44
gradient_avg = 20e6  % [MV/m]

% save(
%     filenameOctaveFieldMap,
%     'R', 'Z', 'Er', 'Ez', 'Btheta', 'Bz',
%     'frequency', 'phase_advance', 'wave_direction', 'cell_length', 'tot_cells', 'gradient_avg');

field1d = load('../../RunningSimulations/RFTrack/YongkeTool_V3/field/field_map_LargeR_Lband.dat');

overallPhaseShift = exp(j*pi*0);

clf();
figure(1);
subplot(5, 1, 1)
plot(field1d.Z, abs(field1d.Ez));
hold on;
plot(Z(:, 1), abs(Ez(:, 1)));
xlabel('Z [m]');
ylabel('abs(Ez) on axis [V/m]');
legend('Old 1D field map (Hermann)', 'New 2D field map (Alexej)');
grid on;
subplot(5, 1, 2)
plot(field1d.Z, real(field1d.Ez*overallPhaseShift));
hold on;
plot(Z(:, 1), real(Ez(:, 1)*overallPhaseShift));
xlabel('Z [m]');
ylabel('Re(Ez) on axis [V/m]');
legend('Old 1D field map (Hermann)', 'New 2D field map (Alexej)');
grid on;
subplot(5, 1, 3)
plot(field1d.Z, imag(field1d.Ez*overallPhaseShift));
hold on;
plot(Z(:, 1), imag(Ez(:, 1)*overallPhaseShift));
xlabel('Z [m]');
ylabel('Im(Ez) on axis [V/m]');
legend('Old 1D field map (Hermann)', 'New 2D field map (Alexej)');
grid on;
subplot(5, 1, 4)
plot(Z(:, 8), real(Btheta(:, 8)*overallPhaseShift), 'r');
xlabel('Z [m]');
ylabel(['Re(Btheta) at r=' num2str(R(1,8)*1e3, 2) 'mm [T]']);
legend('New 2D field map (Alexej)');
grid on;
subplot(5, 1, 5)
plot(Z(:, 8), real(Er(:, 8)*overallPhaseShift), 'r');
xlabel('Z [m]');
ylabel(['Re(Er) at r=' num2str(R(1,8)*1e3, 2) 'mm [T]']);
legend('New 2D field map (Alexej)');
grid on;
figure(2);
subplot(4, 1, 1)
plot(Z(:, 1), abs(Er(:, 1)));
hold on;
plot(Z(:, 8), abs(Er(:, 8)));
plot(Z(:, 16), abs(Er(:, 16)));
xlabel('Z [m]');
ylabel('abs(Er) [V/m]');
legend('On axis', 'Middle r', 'End r');
subplot(4, 1, 2)
plot(Z(:, 1), abs(Btheta(:, 1)));
hold on;
plot(Z(:, 8), abs(Btheta(:, 8)));
plot(Z(:, 16), abs(Btheta(:, 16)));
xlabel('Z [m]');
ylabel('abs(Btheta) [T]');
legend('On axis', 'Middle r', 'End r');
subplot(4, 1, 3)
plot(Z(:, 1), abs(Etheta(:, 1)));
hold on;
plot(Z(:, 8), abs(Etheta(:, 8)));
plot(Z(:, 16), abs(Etheta(:, 16)));
xlabel('Z [m]');
ylabel('abs(Etheta) [V/m]');
legend('On axis', 'Middle r', 'End r');
subplot(4, 1, 4)
plot(Z(:, 1), abs(Br(:, 1)));
hold on;
plot(Z(:, 8), abs(Br(:, 8)));
plot(Z(:, 16), abs(Br(:, 16)));
xlabel('Z [m]');
ylabel('abs(Br) [T]');
legend('On axis', 'Middle r', 'End r');

return
