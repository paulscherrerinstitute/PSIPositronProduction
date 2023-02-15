clear all;

fid = fopen ('pLinacF3_full44cells_Ecomp_YZplane_dy2mm_dz0p1L.fld', 'r');
rawE = textscan(fid, '%f%f%f%f%f%f%f%f%f', 'HeaderLines', 2);
fclose(fid);
rawE = cell2mat(rawE);
size(rawE)
Ex = rawE(:, 4) + j * rawE(:, 5);
Ey = rawE(:, 6) + j * rawE(:, 7);
Ez = rawE(:, 8) + j * rawE(:, 9);

fid = fopen ('pLinacF3_full44cells_Hcomp_YZplane_dy2mm_dz0p1L.fld', 'r');
rawB = textscan(fid, '%f%f%f%f%f%f%f%f%f', 'HeaderLines', 2);
fclose(fid);
rawB = cell2mat(rawB);
size(rawB)
Bx = rawB(:, 4) + j * rawB(:, 5);
By = rawB(:, 6) + j * rawB(:, 7);
Bz = rawB(:, 8) + j * rawB(:, 9);

for axInd = 1:3
    if sum(rawE(:, axInd) ~= rawB(:, axInd)) > 0
        error('Meshes do not correspond!')
    endif
endfor
X = rawE(:, 1);
Y = rawE(:, 2);
Z = rawE(:, 3);

amplitudeScaling = 3.3513e7 / 5.1797e3;
phaseAlignment = exp(pi*j);
matShape = [461, 16];
X = reshape(X, matShape);
Y = reshape(Y, matShape);
Z = reshape(Z, matShape) + 0.11942;
Ex = reshape(Ex, matShape) * amplitudeScaling * phaseAlignment;
Ey = reshape(Ey, matShape) * amplitudeScaling * phaseAlignment;
Ez = reshape(Ez, matShape) * amplitudeScaling * phaseAlignment;
Bx = reshape(Bx, matShape) * amplitudeScaling * phaseAlignment;
By = reshape(By, matShape) * amplitudeScaling * phaseAlignment;
Bz = reshape(Bz, matShape) * amplitudeScaling * phaseAlignment;
size(X)
size(Y)
size(Z)
size(Ex)
size(Ey)
size(Ez)
size(Bx)
size(By)
size(Bz)

frequency = 2e9  % [Hz]
phase_advance = 9./10 * pi  % [rad]
wave_direction = +1
cell_length = 2.99792458e8 / frequency * 9./20.  % [m]
tot_cells = 44
gradient_avg = 20e6  % [MV/m]

% save(
%     'pLinacF3_full44cells_YZplane_dy2mm_dz0p1L.dat',
%     'X', 'Y', 'Z', 'Ex', 'Ey', 'Ez', 'Bx', 'By', 'Bz',
%     'frequency', 'phase_advance', 'wave_direction', 'cell_length', 'tot_cells', 'gradient_avg');

field1d = load('../../RunningSimulations/RFTrack/YongkeTool_V3/field/field_map_LargeR_Lband.dat');

clf();
figure(1);
subplot(3, 1, 1)
plot(field1d.Z, abs(field1d.Ez));
hold on;
plot(Z(:, 1), abs(Ez(:, 1)));
xlabel('Z [m]');
ylabel('Ez [V/m]');
legend('Old 1D field map (Hermann)', 'New 2D field map (Alexej)');
grid on;
subplot(3, 1, 2)
plot(field1d.Z, real(field1d.Ez));
hold on;
plot(Z(:, 1), real(Ez(:, 1)));
xlabel('Z [m]');
ylabel('Re(Ez) [V/m]');
legend('Old 1D field map (Hermann)', 'New 2D field map (Alexej)');
grid on;
subplot(3, 1, 3)
plot(field1d.Z, imag(field1d.Ez));
hold on;
plot(Z(:, 1), imag(Ez(:, 1)));
xlabel('Z [m]');
ylabel('Im(Ez) [V/m]');
legend('Old 1D field map (Hermann)', 'New 2D field map (Alexej)');
grid on;

return
