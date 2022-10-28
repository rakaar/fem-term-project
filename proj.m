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
for f=1:size(force_vecs,2)
    nodal_f = -(1/3)* (500/11)*(f-1);
    force_vecs{1,f}(6,1) = nodal_f; % force applied - 3rd node, 2nd dof
    force_vecs{1,f}(8,1) = nodal_f; % force applied - 3rd node, 2nd dof
    force_vecs{1,f}(14,1) = nodal_f; % force applied - 3rd node, 2nd dof
end

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

%% printing displacement of element
diary 'displacement.txt'
disp('non-zero displacement values - ux, uy of each node')
for i=2:size(displacement_vecs,2)
    fprintf("\n In top of the bar, elemenent num %d \n",i);
    disp(reshape(displacement_vecs{1,i}, 1,16));
end
diary off

%% printing stress at gauss pts
diary 'gauss_stress.txt'
disp('non-zero stress at gauss points - sigma x, sigma y, tau xy')
for i=2:size(displacement_vecs,2)
    fprintf("\n In top of the bar, elemenent num %d \n",i);
    for j=1:4
        disp(gauss_pts_stress_vecs{1,i}{j});
    end
end
diary off

%% printing stress at nodes

diary 'nodal_stress.txt'
disp('non-zero stress at nodes - sigma x, sigma y, tau xy')
for i=2:size(displacement_vecs,2)
    fprintf("\n In top of the bar, elemenent num %d \n",i);
    for j=1:8
        disp(stress_vec_nodes{1,i}{j});
    end
end
diary off


%% plotting ux and uy of all middle node on top side
ux = zeros(12,1);
uy = zeros(12,1);
for i=1:12
    ux(i) = displacement_vecs{1,i}(13);
    uy(i) = displacement_vecs{1,i}(14);
end

figure
    hold on
    plot(ux)
    plot(uy, '--')
    hold off
    title('Ux and Uy of an example node')
    legend('Ux','Uy')
grid
%% plotting average sigma - x, sigma  y, tau - xy of all elements
sigma_x = zeros(12,1);
sigma_y = zeros(12,1);
tau_xy = zeros(12,1);

for i=1:12
    for j=1:8
        sigma_x(i) = sigma_x(i) + stress_vec_nodes{1,i}{j}(1)/8;
        sigma_y(i) = sigma_y(i) + stress_vec_nodes{1,i}{j}(2)/8;
        tau_xy(i) = tau_xy(i) + stress_vec_nodes{1,i}{j}(3)/8;
    end
    
end

figure
    hold on
    plot(sigma_x)
    plot(sigma_y, '--')
    plot(tau_xy, '-o')
    hold off
    title('Average stress in each element - sigma x, sigma y and tau xy')
    legend('sigma x','sigma y', 'tau xy')
grid

