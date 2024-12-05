%%%%% Original file for Problem 1 (incremental)

clear; close all;
%% Mesh generation
elementType='2dQ1'; %define element type

% * ADD: Define boundary element type
elementType1d='1dQ1'; % consistent with elementType
% *

domain=[0 0 4 1]; %[x0 y0 x1 y1]
numElements=[15 10]; %[numElX numElY]
[nodeCoords, IEN, boundaryElementIDs, boundaryNodeLocalID]=...
    meshRect2d(domain,elementType,numElements); %generate the mesh

% * ADD: Create boundary incidence arrays
BIEN=IENtoBIEN(IEN, boundaryElementIDs, boundaryNodeLocalID);
% *
%% Parameters
%elasticity tensor
CMatrix=elasticProperties('youngsModulus',193e6,'poissonsRatio',0.253,'CPlaneStressEng');

%% Matrix and vector assembly 
% ! REMOVE: Remove unused parameters
% numNodes=size(nodeCoords,1); %total number of nodes
% numDoFs=numNodes*2; %number of global degrees of freedom
% !

% Assemble stiffness matrix
numGP=1; %number of Gauss points used in quadrature
K = formStiffnessMatrixEng(nodeCoords, IEN, elementType, numGP, CMatrix);

%define function handle for body force
bodyForce=@(x)(repmat([0; 0],size(x,1))); 
% assemble body force vector
Fb = formBodyForceVector(nodeCoords, IEN, elementType, numGP, bodyForce);

% ! REMOVE: Remove old body force
% % Body force
% Fb = zeros(numDoFs,1);
% !

% * ADD: Boundary conditions
%define function handles for BC
bndDisplacement=@(x)(zeros(size(x))'); % zero displacement
bndTraction0=@(x)(zeros(size(x))'); % zero traction
% Traction vector [1; 0] for the right boundary
bndTraction1=@(x)(repmat([1; 0],size(x,1)));  % non-zero traction
%make cell arrays of function handles - boundary condition cell arrays must
%match the structure of BIEN 
bndTractions={bndTraction0,bndTraction1,bndTraction0,bndTraction0};
bndDisplacements=repmat({bndDisplacement},4,1); 
isDirichlet=[0; 1; 0; 1]; %define which boundary regions have Dirichlet BC

numGP1d=1; %number of Gauss Points for 1d boundary elements
% assemble boundary load vector for Neumann BC
% and evaluate displacements for Dirichlet BC
[u_prescribed, Fs, prescribedDoF, freeDoF]=...
    formBC(nodeCoords,BIEN,elementType1d,numGP1d,...
    bndTractions,bndDisplacements,isDirichlet);
% *

% ! REMOVE: Remove old boundary conditions
% % Boundary conditions
% leftNodes=find(nodeCoords(:,1)==domain(1)); 
% rightNodes=find(nodeCoords(:,1)==domain(3));
% leftXDoF=(leftNodes-1)*2+1;
% leftYDoF=(leftNodes-1)*2+2;
% rightXDoF=(rightNodes-1)*2+1;
% rightYDoF=(rightNodes-1)*2+2;

% prescribedDoF=[leftXDoF; leftYDoF]; %list of prescribed global DoFs
% freeDoF=setdiff(1:numDoFs,prescribedDoF); %list of free DoFs
% u_prescribed=zeros(numDoFs,1); 
% % Surface traction
% Fs = zeros(numDoFs,1); %define global force vector for surface traction
% !
%% Solution 
F=Fb+Fs; %total load vector
%define the free part of load vector
FF=F(freeDoF)-K(freeDoF,prescribedDoF)*u_prescribed(prescribedDoF);
%define the free part of stiffness matrix
KK=K(freeDoF,freeDoF);

% * ADD: Solve linear equations
u=zeros(size(u_prescribed)); %initialise vector of displacements
%solve linear equations
u(freeDoF)=KK\FF;
% *

% ! REMOVE: Remove old solution
% %solve linear equations
% u=zeros(numDoFs,1);
% u(freeDoF)=KK\FF;
% u(prescribedDoF)=u_prescribed(prescribedDoF);
% !
%% Stress recovery 
u2=reshape(u,[2 numel(u)/2])'; %reshape s.t. Ux=u2(:,1), Uy=u2(:,2)
%recover strains at centroids elements
[strain, GPCoords]=recoveryGPEng(u2,nodeCoords,IEN,elementType,1);
%evaluate stresses at the centroids
s2=CMatrix*strain;

%% Visualisation
% * MODIFICATION: Draw the initial undeformed mesh
figure(1);clf;
drawElements(nodeCoords,IEN,elementType); % mesh
hold on;
drawNodes(nodeCoords); % nodes
drawNodes(GPCoords,'x'); % Gauss points
% *

figure(2);clf; 
%draw initial undeformed mesh
drawElements(nodeCoords,IEN,elementType,0,0);
hold on;
factor=1e7; %scaling factor to amplify small deformations
%draw deformed mesh
drawElements(nodeCoords+u2*factor,IEN,elementType,s2(1,:)',.7);
drawNodes(nodeCoords,BIEN{4},{'ks','filled'}); %draw pinned nodes
drawNodes(nodeCoords,BIEN{2},{'ks','filled'}); %draw pinned nodes
title('\sigma_x')