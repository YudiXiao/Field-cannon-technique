Nx = 2180;
Ny = 790;
grid_size = 5e-6;
field = table2array(readtable('A_old.dat'));
%
A_old = zeros(Ny+1,Nx+1);
%
for i=1:1:size(field,1)
    A_old(Ny+1-int16(field(i,2)/grid_size),int16(field(i,1)/grid_size)+1) = field(i,3);
end

field = table2array(readtable('A_new.dat'));
%
A_new = zeros(Ny+1,Nx+1);
%
for i=1:1:size(field,1)
    A_new(Ny+1-int16(field(i,2)/grid_size),int16(field(i,1)/grid_size)+1) = field(i,3);
end

RMS = sqrt((1/((Nx+1)*(Ny+1)))*sum((A_new-A_old).^2,'all'));

% 
field = table2array(readtable('B_x.dat'));
B_x = zeros(Ny+1,Nx+1);
%
for i=1:1:size(field,1)
    B_x(Ny+1-int16(field(i,2)/grid_size),int16(field(i,1)/grid_size)+1) = field(i,3);
end

field = table2array(readtable('B_y.dat'));
B_y = zeros(Ny+1,Nx+1);
%
for i=1:1:size(field,1)
    B_y(Ny+1-int16(field(i,2)/grid_size),int16(field(i,1)/grid_size)+1) = field(i,3);
end

B = sqrt(B_x.^2 + B_y.^2);

figure
imagesc(1:size(B_x,2),1:size(B_x,1),B_x)
title('B_x')
colorbar

figure
imagesc(1:size(B_y,2),1:size(B_y,1),B_y)
title('B_y')
colorbar

figure
imagesc(1:size(B,2),1:size(B,1),B)
title('Magnitude of \it{\bf{B}}')
colorbar

figure
imagesc(1:size(A_new,2),1:size(A_new,1),A_new)
colorbar

%% Function to evaluate y-component
% to determine the optimal length of bump

width_winding = 5.5e-3;
dis_wind = 1.5e-3;

N_wind_x = int16(width_winding/grid_size);
N_dis_wind = int16(dis_wind/grid_size);
N_core_x = 500;
N_core_wind_x = 40;

B_y_winding = zeros(1,N_wind_x + 3);
idx_x = zeros(1,N_wind_x + 3);

for i=1:1:(N_wind_x + 3)
    B_y_winding(1,i) = B_y(Ny+1-N_dis_wind,N_core_x+N_core_wind_x+i-1);
    idx_x(1,i) = double(N_core_x+N_core_wind_x+i-1)*grid_size;
end

%%
B_y_winding_abs = abs(B_y_winding);

% standard deviation of B_y beneath the winding
%load('B_y_abs_1.5mm_dis_0.8mm_bump.mat');
B_y_std = std(B_y_winding_abs);

figure
plot(idx_x,B_y_winding_abs)