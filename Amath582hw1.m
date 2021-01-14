%% HW 1

% looking for unknown signal
% data obtained over 24 hours in half hour increments -> 48 realizations
% source of signal is moving
% need to determine location and path by using the acoustic signature and
% identifying the acoustic admissions. 

%% loading data

clear; close all; clc
datast = load('subdata.mat'); 
data = datast.subdata; % space by time matrix

%% Setup 

L = 10; % spatial (10 and -10 are the max and min x-values)
n = 64; % Fourier modes (# of sinusiods to be added together (# of freqncies) to represent the signal)

x2 = linspace(-L,L,n+1); x = x2(1:n); y=x; z=x;
k = (2*pi/(2*L))*[0:(n/2 - 1) -n/2:-1]; ks = fftshift(k); % freqencies (wave numbers) scaled

[X,Y,Z] = meshgrid(x,y,z);
[Kx,Ky,Kz] = meshgrid(ks,ks,ks);

%% Initial Plotting
for j=1:49
    Un(:,:,:)=reshape(data(:,j),n,n,n); % must there always be the same number of spacial points per dim as n?
    M = max(max(max(abs(Un))));
    close all, isosurface(X,Y,Z,abs(Un)/M,0.7)
    axis([-20 20 -20 20 -20 20]), grid on, drawnow
    pause(1)
end


%% Reshaping into a 3D matric evolving in time
N = 48; % number of realizations
clear Un
for j=1:N
    Un(:,:,:,j)=reshape(data(:,j),n,n,n); 
end
% Each realization of Un is a noisy comlex singal located somewhere in x,y,z space
%% Determining frequency signature through spectrom averaging

Unt = fftn(Un);  % FT of Un

Untave = mean(Unt,4);

% need to find corresping kx ky kz to create filter.

MaxVal = max(max(max(abs(Untave))));

Diffmat = abs(abs(Untave) - MaxVal);

SmallVal = min(min(min(Diffmat)));

% location of peak of avarage signal in frequency domain
KxAve = Kx(abs(Diffmat) == SmallVal);
KyAve = Ky(abs(Diffmat) == SmallVal);
KzAve = Kz(abs(Diffmat) == SmallVal);


%% Filter

% when filter using a unit gaussian the signal should be normalized by max
% M first

% making 3d Gaussian function
sig_x = 1; sig_y = 1; sig_z = 1;
mu_x = max_kx; mu_y = max_ky; mu_z = max_kz;
Gaus3DFilt = @(x,y,z) (1/(4*pi*sig_x*sig_y*sig_z))*exp(-( ...
((x-mu_x).^2)/(2*sig_x^2) + ((y-mu_y).^2)/(2*sig_y^2) + ((z-mu_z).^2)/(2*sig_z^2))); filt_xyz = Gaus3DFilt(Kx,Ky,Kz);
