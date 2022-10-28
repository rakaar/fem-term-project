function [stiffness]=formStiffness2D(GDof,numberElements,...elementNodes,numberNodes,nodeCoordinates,thickness)
% compute stiffness matrix
% for plane stress Q4 elements
stiffness=zeros(GDof);
% 2 by 2 quadrature
[gaussWeights,gaussLocations]=gaussQuadrature('complete');
for e=1:numberElements
%-----------------------------------------------------------------------% %Construction of Matrix of Material Constants
 if e<=96
 E=25e9
 elseif e<=192
 E=50e9
 else
 E=100e9
 end
 poisson=0.3;
 D=E/(1-poisson^2)*[1 poisson 0;poisson 1 0;0 0 (1-poisson)/2];%-----------------------------------------------------------------------indice=elementNodes(e,:); % Connectivity of 'e'th elementelementDof=[ indice indice+numberNodes ]; % Global degree of freedom of 'e'th element [xDof yDof]
ndof=length(indice);
% cycle for Gauss point
for q=1:size(gaussWeights,1)
GaussPoint=gaussLocations(q,:);
xi=GaussPoint(1);
eta=GaussPoint(2);
% shape functions and derivatives
[shapeFunction,naturalDerivatives]=shapeFunctionQ4(xi,eta);% Jacobian matrix, inverse of Jacobian,
 [Jacob,invJacobian,XYderivatives]=...
Jacobian(nodeCoordinates(indice,:),naturalDerivatives);
%-----------------------------------------------------------------------% B matrix
B=zeros(3,2*ndof);
B(1,1:ndof) = XYderivatives(:,1)';
B(2,ndof+1:2*ndof) = XYderivatives(:,2)';
B(3,1:ndof) = XYderivatives(:,2)';
B(3,ndof+1:2*ndof) = XYderivatives(:,1)';
%-----------------------------------------------------------------------% stiffness matrix
stiffness(elementDof,elementDof)= stiffness(elementDof,elementDof)+...B'*D*thickness*B*gaussWeights(q)*det(Jacob);
end
end

end
