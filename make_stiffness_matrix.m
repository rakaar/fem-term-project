function stiffness_matrix = make_stiffness_matrix(element_global_coords)
    % calculates stiffness matrix for all 4 gauss points and returns their
    % sum k = Ka + kb + Kc + Kd
    % input: element_global_coords  8 x 1 cell, each cell is a coordinate in global system

    gauss_pts = [[-1 -1]; [1 -1]; [1 1]; [-1 1]]/(3^0.5);
    gauss_pts_stiffness_matrix = cell(4,1); % Ka,Kb,Kc,Kd
    gauss_pts_g_matrix = cell(4,1); % Ga, Gb, Gc, Gd
    gauss_pts_jacob_matrix = cell(4,1);
    gauss_pts_a_matrix = cell(4,1);
    gauss_pts_b_matrix = cell(4,1);

    % element global coords
    x = zeros(8,1); y = zeros(8,1);
    for i=1:8
        x(i) = element_global_coords{i}(1);
        y(i) = element_global_coords{i}(2);
    end
    
    for g=1:length(gauss_pts)
        xi = gauss_pts(g,1);
        eta = gauss_pts(g,2);
        
        % xi derivatives
        s = zeros(8,1);
        s(1) = 0.25*(1 - eta)*(2*xi + eta);
        s(2) = 0.25*(1 - eta)*(2*xi - eta);
        s(3) = 0.25*(1 + eta)*(2*xi + eta);
        s(4) = 0.25*(1 + eta)*(2*xi - eta);
        s(5) = -xi*(1-eta);
        s(6) = 0.5*(1-(eta^2));
        s(7) = -xi*(1 + eta);
        s(8) = -0.5*(1 - (eta^2));

        % eta derivatives
        t = zeros(8,1);
        t(1) = 0.25*(1 - xi)*(2*eta + xi);
        t(2) = 0.25*(1 + xi)*(xi - 2*eta);
        t(3) = 0.25*(1 + xi)*(2*eta + xi);
        t(4) = 0.25*(1 - xi)*(2*eta - xi);
        t(5) = -0.5*(1 - (xi^2));
        t(6) = -eta*(1 + xi);
        t(7) = 0.5*(1 - (xi^2));
        t(8) = -eta*(1 - xi);

        % g matrix
        g_matrix = [ ...
        [s(1) 0 s(2) 0 s(3) 0 s(4) 0 s(5) 0 s(6) 0 s(7) 0 s(8) 0];...
        [t(1) 0 t(2) 0 t(3) 0 t(4) 0 t(5) 0 t(6) 0 t(7) 0 t(8) 0];...
        [0 s(1) 0 s(2) 0 s(3) 0 s(4) 0 s(5) 0 s(6) 0 s(7) 0 s(8)]; ...
        [0 t(1) 0 t(2) 0 t(3) 0 t(4) 0 t(5) 0 t(6) 0 t(7) 0 t(8)];...
        ];

        gauss_pts_g_matrix{g,1} = g_matrix; 

        % jacobian matrix elements
        J = zeros(2,2);
        for i=1:8
            J(1,1) = J(1,1) + s(i)*x(i);
            J(1,2) = J(1,2) + s(i)*y(i);
            J(2,1) = J(2,1) + t(i)*x(i);
            J(2,2) = J(2,2) + t(i)*y(i);
        end % end of jacobian matrix
        gauss_pts_jacob_matrix{g,1} = J;
        detJ = det(J);

        %  A matrix
        a_matrix = [...
            [J(2,2) -J(1,2) 0 0];...
            [0 0 -J(2,1) J(1,1)];...
            [-J(2,1) J(1,1) J(2,2) -J(1,2)];...
        ];
        gauss_pts_a_matrix{g,1} = a_matrix;

        % B matrix
        b_matrix = a_matrix*g_matrix;
        gauss_pts_b_matrix{g,1} = b_matrix;

        % D matrix 
        % matrix ????
        
        % K matrix
        % B.' D B
    end % end of for each gauss pt
end