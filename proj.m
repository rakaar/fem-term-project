clear all
%% physical information
length_beam = 96; % in mm
width_beam = 16; % in mm

Young_modulus_1 = 100e9; % in Pa
Young_modulus_2 = 50e9; % in Pa
Young_modulus_3 = 25e9; % in Pa

F = 500; % in N, force at end of beam

%% Gauss quadrature
xi = [-(3/5)^(0.5) 0 (3/5)^(0.5)]; % points
wi = [5/9 8/9 5/9]; % weights

%% fem nodes - 8 noded element
node_size = 4; % in mm
element_vertices = 8;

x = 0:node_size:length_beam;
y = 0:node_size:3*width_beam;

[X,Y] = meshgrid(x,y);

nodal_coords = cell(size(X,1),size(X,2)); % beam top edge is on X-axis
for i=1:size(X,1)
    for j=1:size(X,2)
        nodal_coords{i,j} = [X(i,j) Y(i,j)];
    end
end

element_coords = cell(size(nodal_coords,1) - 2, size(nodal_coords,2) - 2);
for r=1:size(nodal_coords,1) - 2
    for c=1:size(nodal_coords,2) - 2
        element_coords{r,c} = [nodal_coords(r+2,c); nodal_coords(r+2,c+2); nodal_coords(r,c+2); nodal_coords(r,c); nodal_coords(r+2,c+1); nodal_coords(r+1,c+2); nodal_coords(r,c+1); nodal_coords(r+1,c)];
    end
end



stiffness_matrices = cell(size(nodal_coords,1) - 2, size(nodal_coords,2) - 2);
for r=1:size(nodal_coords,1) - 2
    for c=1:size(nodal_coords,2) - 2
        stiffness_matrices{r,c} = make_stiffness_matrix(element_coords{r,c});
    end

end
