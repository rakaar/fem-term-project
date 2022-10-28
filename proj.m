clear all
%% physical information
length_beam = 96; % in mm
width_beam = 16; % in mm

Young_modulus_1 = 100e9; % in Pa
Young_modulus_2 = 50e9; % in Pa
Young_modulus_3 = 25e9; % in Pa
poisson_ratio = 0.3;

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

element_coords = cell((size(nodal_coords,1) -1)/2, (size(nodal_coords,2) - 1)/2);
for r=1:2:size(nodal_coords,1) - 2
    for c=1:2:size(nodal_coords,2) - 2
        element_coords{(r-1)/2 + 1,(c-1)/2 + 1} = [nodal_coords(r+2,c); nodal_coords(r+2,c+2); nodal_coords(r,c+2); nodal_coords(r,c); nodal_coords(r+2,c+1); nodal_coords(r+1,c+2); nodal_coords(r,c+1); nodal_coords(r+1,c)];
    end
end



stiffness_matrices = cell(size(element_coords,1), size(element_coords,2));
gauss_pts_b_matrices = cell(size(element_coords,1), size(element_coords,2));
for r=1:size(element_coords,1)
    if r >= 1 && r <= 2
        young_modulus = Young_modulus_1;
    elseif r >= 3 && r <= 4
        young_modulus = Young_modulus_2;
    elseif r >= 5
        young_modulus = Young_modulus_3;
    end
    
    for c=1:size(element_coords,2)
[stiffness_matrices{r,c},gauss_pts_b_matrices{r,c}] = make_stiffness_matrix(element_coords{r,c}, young_modulus, poisson_ratio);
    end
end

%% force vectors for all elements
force_vecs = cell(size(element_coords,1), size(element_coords,2)); % each is 16x1 - [Fu1 Fv1 Fu2 Fv2 ... Fu8 Fv8]
for r=1:size(element_coords,1)
    for c=1:size(element_coords,2)
        force_vecs{r,c} = zeros(16,1);
    end
end

% except element at the right corner
force_vecs{1,size(force_vecs,2)}(6,1) = -F/3; % force applied - 3rd node, 2nd dof
force_vecs{1,size(force_vecs,2)}(8,1) = -F/3; % force applied - 3rd node, 2nd dof
force_vecs{1,size(force_vecs,2)}(14,1) = -F/3; % force applied - 3rd node, 2nd dof
%% displacement F = KQ
displacement_vecs = cell(size(element_coords,1), size(element_coords,2)); % each is 16x1 - [u1 v1 u2 v2 ... u8 v8]
for r=1:size(element_coords,1)
   
    for c=1:size(element_coords,2)
        stiffness_matrices{r,c} = make_stiffness_matrix(element_coords{r,c}, young_modulus, poisson_ratio);
        displacement_vecs{r,c} = pinv(stiffness_matrices{r,c})*force_vecs{r,c};
    end
end

%% calculate stress at gauss pts
gauss_pts_stress_vecs = cell(size(element_coords,1), size(element_coords,2));

for r=1:size(element_coords,1)
    for c=1:size(element_coords,2)
             if r >= 1 && r <= 2
                young_modulus = Young_modulus_1;
            elseif r >= 3 && r <= 4
                young_modulus = Young_modulus_2;
            elseif r >= 5
                young_modulus = Young_modulus_3;
             end

       d_matrix = (young_modulus/(1 - poisson_ratio^2))* [...
            [1 poisson_ratio 0];...
            [poisson_ratio 1 0];...
            [0 0 (1-poisson_ratio)/2]
        ];
          stress_cell = cell(4,1);
          for g=1:4
              stress_cell{g,1} = d_matrix*gauss_pts_b_matrices{r,c}{g}*displacement_vecs{r,c};
          end

          gauss_pts_stress_vecs{r,c} = stress_cell;
          
    end
end

%% calculate stress at nodes
stress_vec_nodes = cell(size(element_coords,1), size(element_coords,2));
for r=1:size(element_coords,1)
    for c=1:size(element_coords,2)
        if r >= 1 && r <= 2
                young_modulus = Young_modulus_1;
            elseif r >= 3 && r <= 4
                young_modulus = Young_modulus_2;
            elseif r >= 5
                young_modulus = Young_modulus_3;
             end

       d_matrix = (young_modulus/(1 - poisson_ratio^2))* [...
            [1 poisson_ratio 0];...
            [poisson_ratio 1 0];...
            [0 0 (1-poisson_ratio)/2]
        ];
  
            stress_vec_nodes{r,c} = calc_stress_at_nodes(d_matrix,displacement_vecs{r,c}, element_coords{r,c});
    end
end

