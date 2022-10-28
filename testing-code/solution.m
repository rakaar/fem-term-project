Function[displacements]=solution(GDof,prescribedDof,...
 Stiffness, force)
%freeDof= Active Dof
freeDof=setdiff((1:GDof)',prescribedDof);
displacements(freeDof)=inv(stiffness(freeDof,freeDof))*(force(freeDof));