% MATLAB codes for Finite Element Analysis
% Plane Stress problem
%---------------------------------------------------------------- % clear memory
clear all;close all; clc ;
colordef white;clf
format short;
%----------------------------------------------------------------% load
P = -500;
%----------------------------------------------------------------%Mesh generation
Lx=96; Ly=48;
numberElementsX=24;
numberElementsY=12;
numberElements=numberElementsX*numberElementsY;
[nodeCoordinates, elementNodes] = rectangularMesh(Lx,Ly,numberElementsX,numberElementsY);
xx=nodeCoordinates(:,1);
yy=nodeCoordinates(:,2);
numberNodes = size(xx,1);
GDof = 2*numberNodes; % GDof: global number of degrees of freedomPlotMesh(nodeCoordinates, elementNodes);
axis equal
%----------------------------------------------------------------% calculation of the system stiffness matrix
stiffness = formStiffness2D(GDof,numberElements,elementNodes,numberNodes,nodeCoordinates,D,1);
%----------------------------------------------------------------% boundary conditions
fixedNodeX = find(nodeCoordinates(:,1)==0); % fixed in XXfixedNodeY = find(nodeCoordinates(:,1)==0); % fixed in YYprescribedDof = [fixedNodeX; fixedNodeY+numberNodes];
%----------------------------------------------------------------% force vector (Concentrated Load at xx=Lx, yy=Ly)
force=zeros(GDof,1);
rightend=find(nodeCoordinates(:,1)==Lx);
force(rightend+numberNodes)= P;
%----------------------------------------------------------------% Solution
displacements = solution(GDof,prescribedDof,stiffness,force); %----------------------------------------------------------------%Separation of X and Y Displacemets
jj = 1:GDof/2;
UX = displacements(1:numberNodes);
UY = displacements(numberNodes+1:GDof);
%----------------------------------------------------------------% Stresses and Strains(at Gauss-points)
[stress,strain] = stresses2D(GDof,numberElements,...
elementNodes,numberNodes,nodeCoordinates,displacements,...
D,'Q4','complete'); 
% 
% %-----------------------------------------------------------------------% Calculation of Stresses and Strains at Node points
% %----------------------------------------------------------------------- % STRAIN AND stress extrapolation
% stressExtr = zeros(numberElements,4,3);
% strainExtr = zeros(numberElements,4,3);
% for e = 1:numberElements
% for i = 1:3
%  stressExtr(e,:,i) =[1+0.5*sqrt(3) -0.5 1-0.5*sqrt(3) -0.5;  -0.5 1+0.5*sqrt(3) -0.5 1-0.5*sqrt(3);  1-0.5*sqrt(3) -0.5 1+0.5*sqrt(3) -0.5;  -0.5 1-0.5*sqrt(3) -0.5
% 1+0.5*sqrt(3)]*...
%  [stress(e,1,i);stress(e,2,i);stress(e,3,i);stress(e,4,i)]; 
%  strainExtr(e,:,i) =[1+0.5*sqrt(3) -0.5 1-0.5*sqrt(3) -0.5;  -0.5 1+0.5*sqrt(3) -0.5 1-0.5*sqrt(3);  1-0.5*sqrt(3) -0.5 1+0.5*sqrt(3) -0.5;  -0.5 1-0.5*sqrt(3) -0.5
% 1+0.5*sqrt(3)]*...
%  [strain(e,1,i);strain(e,2,i);strain(e,3,i);strain(e,4,i)]; end
% end
% %-----------------------------------------------------------------------% stress and strain averaging at nodes
% stressAvg = zeros(numberNodes,3);
% strainAvg = zeros(numberNodes,3);
% for i = 1:3
% currentStress = stressExtr(:,:,i);
% currentStrain = strainExtr(:,:,i);
% for n = 1:numberNodes
% idx = find(n==elementNodes);
% stressAvg(n,i) = sum(currentStress(idx))/...
% length(currentStress(idx));
% strainAvg(n,i) = sum(currentStrain(idx))/...
% length(currentStrain(idx));
% end
% end
% %----------------------------------------------------------------% POST- PROCESSING
% %----------------------------------------------------------------%Printing values of Field variables at nodes in text format diary Displacements.txt
% table(jj',UX',UY')
% diary off
% diary Strains.txt
% table(jj',strainAvg(:,1).*10e6,strainAvg(:,2).*10e6,...
% strainAvg(:,3).*10e6)
% diary off
% diary Stress.txt
% table(jj',stressAvg(:,1),stressAvg(:,2),stressAvg(:,3))
% diary off
% %Displacment Plots on Mesh of the beam
% PlotFieldonMesh(nodeCoordinates, elementNodes,UY');
% axis equal; xlabel('X(mm)'); ylabel ('Y(mm)');
% title('Variation of Nodal Displacement(in mm) along vertical direction'); 
% 
% PlotFieldonMesh(nodeCoordinates, elementNodes,UX');
% axis equal; xlabel('X(mm)'); ylabel ('Y(mm)');
% title('Variation of Nodal Displacement(in mm) along horizontal direction');
% 
% %Strain Plots on Mesh of the beam
% PlotFieldonMesh(nodeCoordinates, elementNodes,strainAvg(:,1)); axis equal; xlabel('X(mm)'); ylabel ('Y(mm)');
% title('Variation of Normal Strain component along X direction'); 
% PlotFieldonMesh(nodeCoordinates, elementNodes,strainAvg(:,2)); axis equal; xlabel('X(mm)'); ylabel ('Y(mm)');
% title('Variation of Normal Strain component along Y direction'); 
% %Stress Plots on Mesh of the beam
% PlotFieldonMesh(nodeCoordinates, elementNodes,stressAvg(:,1)); axis equal; xlabel('X(mm)'); ylabel ('Y(mm)');
% title('Variation of Normal Stress component(in MPa) along X direction'); 
% PlotFieldonMesh(nodeCoordinates, elementNodes,stressAvg(:,2)); axis equal; xlabel('X(mm)'); ylabel ('Y(mm)');
% title('Variation of Normal Stress component (in MPa) along Y direction'); 