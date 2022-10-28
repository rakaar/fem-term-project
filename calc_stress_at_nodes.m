function stress_vecs = calc_stress_at_nodes(d_matrix,displacement_vec,global_element_coords)
    % returns stress vec at each node of an element
    % given an element's nodes' global coords, D matrix(Element Constant Matrix), Displacement vector      
    stress_vecs = cell(8,1); % eight 3x1 vecs
    xi_vec = zeros(8,1); eta_vec = zeros(8,1);
    xi_vec(1) = -1; xi_vec(2) = 1; xi_vec(3) = 1; xi_vec(4) = -1; xi_vec(5) = 0; xi_vec(6) = 1; xi_vec(7) = -1; xi_vec(8) = 0;
    eta_vec(1) = -1; eta_vec(2) = -1; eta_vec(3) = 1; eta_vec(4) = 1; eta_vec(5) = -1; eta_vec(6) = 0; eta_vec(7) = 1; eta_vec(8) = 0;
     x = zeros(8,1); y = zeros(8,1);
   
    for i=1:8
        x(i) = global_element_coords{i}(1);
        y(i) = global_element_coords{i}(2);
    end

    for n=1:8
        xi = xi_vec(n);
        eta = eta_vec(n);

        % build g matrix
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

        % build 'a' matrix
        % jacobian matrix elements
        J = zeros(2,2);
        for i=1:8
            J(1,1) = J(1,1) + s(i)*x(i);
            J(1,2) = J(1,2) + s(i)*y(i);
            J(2,1) = J(2,1) + t(i)*x(i);
            J(2,2) = J(2,2) + t(i)*y(i);
        end % end of jacobian matrix
        

        %  A matrix
        a_matrix = [...
            [J(2,2) -J(1,2) 0 0];...
            [0 0 -J(2,1) J(1,1)];...
            [-J(2,1) J(1,1) J(2,2) -J(1,2)];...
        ];

        stress_vecs{n,1} = d_matrix*a_matrix*g_matrix*displacement_vec;
    end

end

