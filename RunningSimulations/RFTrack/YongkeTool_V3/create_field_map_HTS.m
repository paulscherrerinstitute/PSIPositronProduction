
infilename = 'field/field_map_HTS_5coils_Apr2022_Original.txt';

A = dlmread(infilename);

% Remove first row

A = A(2:end,1:4);

% Turn	[ R/m Z/m Bz/T Br/T ]
% into	[ R/mm Z/mm Node (not used) Br/T Bz/T ]

B(:,1) = A(:,1) * 1e3;
B(:,2) = A(:,2) * 1e3;
B(:,4) = A(:,4);
B(:,5) = A(:,3);
A = B;

% Remove first & last data

zmin = min(A(:,2)); %% mm
zmax = max(A(:,2)); %% mm

M = A(:,2) != zmin & A(:,2) != zmax;
A = A(M,:);

% Get range

zmin = min(A(:,2)); %% mm
zmax = max(A(:,2)); %% mm
rmax = max(A(:,1)); %% mm

% Steps

dr = 0.1; % mm, grid size
dz = 0.1; % mm, grid size

% AMD parameters

amd.R = rmax; %% [mm]
amd.Zmin = zmin; %% [mm]
amd.Zmax = zmax; %% [mm]
amd.L_total = amd.Zmax - amd.Zmin; %% [mm]

printf("INFO:: Field step [dr,dz]: [%f, %f] [mm]\n",dr,dz);
printf("INFO:: Field R range: [0, %f] [mm]\n",amd.R);
printf("INFO:: Field Z range: [%f, %f] [mm]\n",amd.Zmin,amd.Zmax);
printf("INFO:: Total field length: %f [mm]\n",amd.L_total);

nr = floor(amd.R/dr); % grid points
nz = floor(amd.L_total/dz); % grid points

Zaxis = linspace(amd.Zmin, amd.Zmax, nz+1); % mm (1D variable)
Raxis = linspace(0, amd.R, nr+1); % mm (1D variable)

[Z,R] = meshgrid(Zaxis, Raxis); % (2D var.)

Br = griddata(A(:,2), A(:,1), A(:,4), Z, R, 'linear'); % (2D)
Bz = griddata(A(:,2), A(:,1), A(:,5), Z, R, 'linear'); % (2D)

outfilename = 'field/field_map_HTS_5coils_Apr2022_Processed.dat';
save ("-text", outfilename, "Z","R","Br","Bz"); % (2D)

