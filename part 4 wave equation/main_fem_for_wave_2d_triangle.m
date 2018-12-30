clear
time_step = 0.01;
k = 1;
xl = 0; xr = 1; yl = 0; yr = 1;
u_D=@(x) x(:,1)*0;
N = 2^6;
[c4n,n4e,ind4e,n4db] = mesh_fem_2d_triangle(xl,xr,yl,yr,N,N,k);
[M_R, Srr_R, Srs_R, Ssr_R, Sss_R, Dr_R, Ds_R] = get_matrices_2d_triangle(k);
number_of_nodes = size(c4n,1);
u0 = zeros(number_of_nodes,1);
u1 = zeros(number_of_nodes,1);
for time = 0:time_step:1
    f = @(x) time*2*pi^2*sin(pi*x(:,1)).*sin(pi*x(:,2));
    temp_u = fem_for_wave_2d_triangle(c4n,n4e,n4db,ind4e,M_R,Srr_R,Srs_R, Ssr_R,Sss_R,f,u_D, time_step, u1, u0);
    u0 = u1;
    u1 = temp_u;
end
