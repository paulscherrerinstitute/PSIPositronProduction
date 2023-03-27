clear all;

fid = fopen ('field_map_chicane_all.fld', 'r');
rawB = textscan(fid, '%f%f%f%f%f%f', 'HeaderLines', 2);
fclose(fid);
rawB = cell2mat(rawB);
size(rawB)

rawB = sortrows(rawB, [3, 2, 1]);

X = rawB(:, 1);
Y = rawB(:, 2);
Z = rawB(:, 3);
Bx = rawB(:, 4);
By = rawB(:, 5);
Bz = rawB(:, 6);

figure(1);
clf();
subplot(3, 1, 1);
plot(X);
xlabel('Matrix index');
ylabel('X [m]');
grid on;
subplot(3, 1, 2);
plot(Y);
xlabel('Matrix index');
ylabel('Y [m]');
grid on;
subplot(3, 1, 3);
plot(Z);
xlabel('Matrix index');
ylabel('Z [m]');
grid on;

matShape = [65, 9, 601];
X = reshape(X, matShape);
Y = reshape(Y, matShape);
Z = reshape(Z, matShape);
Bx = reshape(Bx, matShape);
By = reshape(By, matShape);
Bz = reshape(Bz, matShape);
size(X)
size(Y)
size(Z)
size(Bx)
size(By)
size(Bz)

field_peak = 0.197  % [T]

% save('field_map_chicane_all.dat', 'X', 'Y', 'Z', 'Bx', 'By', 'Bz', 'field_peak');

figure(2);
clf();
plot(squeeze(Z(33, 5, :)), squeeze(By(33, 5, :)), 'k-');
hold on;
plot(squeeze(Z(33, 1, :)), squeeze(By(33, 1, :)), '--');
plot(squeeze(Z(33, 9, :)), squeeze(By(33, 9, :)), ':');
plot(squeeze(Z(1, 5, :)), squeeze(By(1, 5, :)), '--');
plot(squeeze(Z(65, 5, :)), squeeze(By(65, 5, :)), ':');
xlabel('Z [m]');
ylabel('By [T]');
legend('On-axis', '@ y = -20 mm', '@ y = +20 mm', '@ x = -160 mm', '@ x = +160 mm');
grid on;

return
