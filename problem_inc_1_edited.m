%%%%% Original file for Problem 1 (incremental)
clear; close all;
%% Mesh generation
elementType='2dQ1'; %define element type
domain=[0 0 4 1]; %[x0 y0 x1 y1]
numElements=[15 10]; %[numElX numElY]
[nodeCoords, IEN, boundaryElementIDs, boundaryNodeLocalID]=...
    meshRect2d(domain,elementType,numElements); %generate the mesh

%% Parameters
%elasticity tensor
CMatrix=elasticProperties('youngsModulus',193e6,'poissonsRatio',0.253,'CPlaneStressEng');

%% Matrix and vector assembly 
numNodes=size(nodeCoords,1); %total number of nodes
numDoFs=numNodes*2; %number of global degrees of freedom
% Assemble stiffness matrix
numGP=4; %number of Gauss points used in quadrature
K = formStiffnessMatrixEng(nodeCoords, IEN, elementType, numGP, CMatrix);
% Body force
Fb = zeros(numDoFs,1);

% Boundary conditions
leftNodes=find(nodeCoords(:,1)==domain(1)); % Nodes on the left edge (x = 0)
rightNodes=find(nodeCoords(:,1)==domain(3)); % Nodes on the right edge (x = L)
leftXDoF=(leftNodes-1)*2+1; % Left edge x-direction DoFs
leftYDoF=(leftNodes-1)*2+2; % Left edge y-direction DoFs
rightXDoF=(rightNodes-1)*2+1; % Right edge x-direction DoFs
rightYDoF=(rightNodes-1)*2+2; % Right edge y-direction DoFs

% Begin with DoF that are fixed/zero 
prescribedDoF=[leftXDoF; leftYDoF]; %list of prescribed global DoFs

% Set the left boundary & right boundary y displacement to zero
K(prescribedDoF,:) = 0;
K(:, prescribedDoF) = 0; % Zero out rows and columns
K(prescribedDoF, prescribedDoF) = eye(numel(prescribedDoF)); % Put ones in diagonal

% Add in DoF that are non-zero 
% Combine all prescribed DoFs
prescribedDoF = [leftXDoF; leftYDoF; rightXDoF; rightYDoF]; % Fix both x and y on left and right edges

% Initialise prescribed displacement vector 
u_prescribed=zeros(numDoFs,1); 

% Apply a stretching displacement at the right boundary (x direction only)
stretchDisplacement = 1e-7; % Example stretch displacement
u_prescribed(rightXDoF) = stretchDisplacement; % Apply the stretch displacement to the right

% Penalty method to enforce boundary conditions 
% penalty = 1e10; % large number 
% K(prescribedDoF,prescribedDoF) = K(prescribedDoF, prescribedDoF) + penalty * eye(numel(prescribedDoF));
% Reduced conditioning number but did not fix the output stresses being a
% bit off

freeDoF=setdiff(1:numDoFs,prescribedDoF); %list of free DoFs

% Body Force
% Apply an upward force at right boundary 
upwardLoad = 1e-2; % Example upward force
% upwardLoad = 0; % Example upward force
Fb = zeros(numDoFs,1);
% Fb = Fb + upwardLoad; % Apply force to all nodes

% Surface traction
Fs = zeros(numDoFs,1) + upwardLoad; %define global force vector for surface traction

disp(Fs)

% Check conditioning number of the stiffness matrix 
condNumber = cond(K);

%% Solution 
F=Fb+Fs; %total load vector
%define the free part of load vector
FF=F(freeDoF)-K(freeDoF,prescribedDoF)*u_prescribed(prescribedDoF);
%define the free part of stiffness matrix
KK=K(freeDoF,freeDoF);
%solve linear equations
u=zeros(numDoFs,1);
u(freeDoF)=KK\FF;
u(prescribedDoF)=u_prescribed(prescribedDoF);

%% Stress recovery 
u2=reshape(u,[2 numel(u)/2])'; %reshape s.t. Ux=u2(:,1), Uy=u2(:,2)
%recover strains at centroids elements
[strain, GPCoords]=recoveryGPEng(u2,nodeCoords,IEN,elementType,1);
%evaluate stresses at the centroids
s2=CMatrix*strain;

%% Visualisation 
figure(1);clf; 
%draw initial undeformed mesh
drawElements(nodeCoords,IEN,elementType,0,0);
hold on;
factor=1e7; %scaling factor to amplify small deformations
%draw deformed mesh
drawElements(nodeCoords+u2*factor,IEN,elementType,s2(1,:)',.7);
title('\sigma_x')

% Tests
% disp(u(leftXDoF)); % Should be 0
% disp(u(leftYDoF)); % Should be 0
% disp(u(rightXDoF)); % Should match stretchDisplacement
% disp(u(rightYDoF)); % Should be 0