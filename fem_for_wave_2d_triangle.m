function [u, A, M, b, fns] = fem_for_wave_2d_triangle(c4n,n4e,n4db, ...
    ind4e,M_R,Srr_R,Srs_R,Ssr_R,Sss_R,f,u_D, time_step, pre_u1, pre_u0)
%fem_for_heat_2d_triangle   FEM solver for wave equation problem in 2D with triangular elements
%   Parameters:
%     - c4n : coordinates for nodes.
%     - n4e : nodes for elements.
%     - n4db : nodes for Dirichlet boundary.
%     - ind4e : indices for elements
%     - M_R : Mass matrix on the reference triangle
%     - Srr_R : Stiffness matrix on the reference triangle
%     - Srs_R : Stiffness matrix on the reference triangle
%     - Ssr_R : Stiffness matrix on the reference triangle
%     - Sss_R : Stiffness matrix on the reference triangle
%     - f : RHS in the Heat Equation
%     - u_D : Dirichlet boundary condition for the solution u
%     - time_step: discretisaation for time
%     - pre_u1: u at the previous time step
%     - pre_u0: u at the previous previous time step
%
%   Returns:
%     - u : numerical solution
%     - A : Global stiffness matrix
%     - M : Global mass matrix
%     - b : Global right-hand side
%     - fns : free nodes

number_of_nodes = size(c4n,1);
A = sparse(number_of_nodes,number_of_nodes);
M = sparse(number_of_nodes,number_of_nodes);

b = zeros(number_of_nodes,1);
u = b;
for j=1:size(n4e,1)
    xr = (c4n(n4e(j,1),1)-c4n(n4e(j,3),1))/2; 
    yr = (c4n(n4e(j,1),2)-c4n(n4e(j,3),2))/2;
    xs = (c4n(n4e(j,2),1)-c4n(n4e(j,3),1))/2; 
    ys = (c4n(n4e(j,2),2)-c4n(n4e(j,3),2))/2;
    J = xr*ys-xs*yr;
    rx=ys/J; ry=-xs/J; sx=-yr/J; sy=xr/J;

    A(ind4e(j,:),ind4e(j,:)) = A(ind4e(j,:),ind4e(j,:)) ...
        + J*((rx^2+ry^2)*Srr_R + (rx*sx+ry*sy)*(Srs_R+Ssr_R) + (sx^2+sy^2)*Sss_R);
    M(ind4e(j,:),ind4e(j,:)) = M(ind4e(j,:),ind4e(j,:)) ...
        + J*M_R;

    b(ind4e(j,:)) = b(ind4e(j,:)) + J*M_R*f(c4n(ind4e(j,:),:));
end
fns = setdiff(1:number_of_nodes, n4db);
u(n4db) = u_D(c4n(n4db,:));
u(fns) = M(fns,fns)\(time_step * time_step * b(fns) - time_step * time_step * A(fns, fns) * pre_u1(fns) + 2.0 * M(fns,fns)*pre_u1(fns) - M(fns,fns)*pre_u0(fns));
end
