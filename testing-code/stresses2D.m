function [stress,strain] = stresses2D(GDof,numberElements, ...elementNodes,numberNodes,nodeCoordinates, ...
displacements,elemType,quadType)
% quadrature according to quadType
[gaussWeights,gaussLocations] = gaussQuadrature(quadType);% stresses at nodes
strain = zeros(numberElements,size(gaussLocations,1),3);
stress = zeros(numberElements,size(gaussLocations,1),3);
% stressPoints = [-1 -1;1 -1;1 1;-1 1];
for e = 1:numberElements
%-----------------------------------------------------------------------% %Construction of Matrix of Material Constants
 if e<=96;
 E=25e9;
 elseif e<=192;
 E=50e9;
 else
 E=100e9;
 end
 poisson=0.3;
 D=E/(1-poisson^2)*[1 poisson 0;poisson 1 0;0 0 (1-poisson)/2];%----------------------------------------------------------------------- indice = elementNodes(e,:);
 elementDof = [ indice indice+numberNodes ];
 nn = length(indice);
 for q = 1:size(gaussWeights,1)
 pt = gaussLocations(q,:);
 wt = gaussWeights(q);
 xi = pt(1);
 eta = pt(2);
% shape functions and derivatives
[shapeFunction,naturalDerivatives] = ...
shapeFunctionQ4(xi,eta);
% Jacobian matrix, inverse of Jacobian,
% derivatives w.r.t. x,y
[Jacob,invJacobian,XYderivatives] = ...
Jacobian(nodeCoordinates(indice,:),naturalDerivatives);
% B matrix
B = zeros(3,2*nn);
B(1,1:nn) = XYderivatives(:,1)';
B(2,nn+1:2*nn) = XYderivatives(:,2)';
B(3,1:nn) = XYderivatives(:,2)';
B(3,nn+1:2*nn) = XYderivatives(:,1)';
% element deformation
STRAIN = B*displacements(elementDof)';
strain(e,q,:)=STRAIN;
stress(e,q,:) = D*STRAIN;
end
end
end