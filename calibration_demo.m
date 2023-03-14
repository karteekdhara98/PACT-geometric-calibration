%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code is part of the research paper: "A method for the geometric 
% calibration of ultrasound transducer arrays with arbitrary geometries".
% 
% The paper is authored by: Karteekeya Sastry, Yang Zhang, Peng Hu, 
% Yilin Luo, Xin Tong, Shuai Na, Lihong V. Wang. The paper is currently 
% in preparation.
%
% Author Affiliation: Caltech Optical Imaging Laboratory, Department of 
% Electrical Engineering, Department of Medical Engineering, California 
% Institute of Technology, 1200 East California Boulevard, Pasadena, 
% CA 91125, USA.
%
% For questions or comments about the code, please contact Karteekeya
% Sastry at sdharave \at caltech.edu .
%
% This code is licensed under the MIT license. 
% See the LICENSE file for details.
%
% Author name: Karteekeya Sastry
% Caltech Optical Imaging Laboratory
% Date: 14 March, 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description of the code:
% This is a demo of the geometric calibration procedure described in the 
% above-mentioned paper. It estimates the coordinates of a simulated 
% arc-shaped array using time-of-arrival (ToA) data from 27 point sources
% in a 3 x 3 x 3 arrangement. The code is divided into 5 steps. 
% 1. Initialize the experiment parameters and load the ToA data.
% 2. Initialize the point source locations. 
% 3. Populate the A matrix and b vector as shown in Section 2 of the paper. 
% 4. Estimate the transducer locations.
% 5. Visualize the results.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars
close all
clc

tic;

%% 1. Initialize experiment parameters and load ToA data
fprintf('1. Initialize parameters: ');
% Speed of sound (SoS)
c = 1500; % m/s, to be obtained from temperature measurements

% Load the ToAs
load('ToAs.mat','ToAs'); % ToAs must be of the size N_pt x N_tx, where N_pt 
% is the number of point sources and N_tx is the number of transducers.
dists = c*ToAs; % distances calculated from the  ToAs and the SoS

% Point source parameters
dx = 1e-2; % m, pitch of point sources in the x direction
dy = 1e-2; % m, pitch of point sources in the y direction
dz = 1e-2; % m, pitch of point sources in the z direction

Nx = 3; % Number of point sources in the x direction
Ny = 3; % Number of point sources in the y direction
Nz = 3; % Number of point sources in the z direction

N_pt = Nx*Ny*Nz; % Total number of point sources
fprintf('Done.\n')

%% 2. Initialize the point source locations
% Create Nx x Ny x Nz grid of points
fprintf('2. Initialize point source locations: ')
x = linspace(0,(Nx-1)*dx,Nx);
y = linspace(0,(Ny-1)*dy,Ny);
z = linspace(0,(Nz-1)*dz,Nz);
[X,Y,Z] = meshgrid(x,y,z);

points = [X(:), Y(:), Z(:)]; % N_pt x 3 array of point source locations
fprintf('Done.\n')

%% 3. Populate A matrix and b vectors (defined in Section 2 of the paper)
fprintf('3. Populate A matrix and b vectors: ')
A = zeros((N_pt*(N_pt-1))/2,3);
b = zeros(size(dists,1),(N_pt*(N_pt-1))/2); % b vector for each transducer
k = 1;

% Loop through all N_pt*(N_pt-1)/2 combinations of the point sources 
for i = 1:N_pt
    point1 = points(i,:); 
    for j = i+1:N_pt
        point2 = points(j,:);
        A(k,:) = point1 - point2; % Each row of A is simply the difference 
        % between the point source coordinates
        tmp = norm(point1)^2-norm(point2)^2;
        b(:,k) = (dists(:,j).^2-dists(:,i).^2+tmp)/2;
        k = k+1;
    end
end
fprintf('Done.\n')

%% 4. Solve for transducer positions
fprintf('4. Solve for transducer positions: ')
tx_x = zeros(size(dists,1),3); % Variable to hold the estimated transducer locations

% For larger transducer array sizes, this loop can be parallelized
for m = 1:size(dists,1)
    tx_x(m,:) = (A\b(m,:)')';
end
fprintf('Done.\n')
toc;

%% 5. Visualize the results
figure; 
scatter3(points(:,1),points(:,2),points(:,3),'ro'); hold on;
scatter3(tx_x(:,1),tx_x(:,2),tx_x(:,3),'b.');
legend('Point sources','Estimated transducer locations');
title('Calibrated transducer locations');
xlabel('x'); ylabel('y'); zlabel('z');
