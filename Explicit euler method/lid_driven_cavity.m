%Initializing parameters
Re = 27;
Nx = 51; Ny = 21;
Lx = 3; Ly = 1;

x = linspace(0,Lx,Nx); dx = x(2) -x(1);
y = linspace(0,Ly,Ny); dy = y(2) - y(1);

[x2,y2] = meshgrid(x,y);

dt = 0.005; tf = 2; Nmax = tf/dt;

%Generating A matrix
[A] = assembleA(Nx,Ny,dx,dy);

%Initializing vort matrix
vort = zeros(Nx-2,Ny-2);
t = 0;

%Solving for streamfunction
for i = 1:Nmax
    stmfunc = solve_Poisson(vort,A,Nx,Ny);
    vort = advanceVort(stmfunc,vort,Nx,Ny,dx,dy,dt,Re,t);
    t = t + dt
end

%Displaying final uv
stmfunc = [0, zeros(1,Ny-2), 0;
            zeros(Nx-2,1), stmfunc, zeros(Nx-2,1);
            0, zeros(1,Ny-2), 0];
u = (stmfunc(2:end-1, 3:end) - stmfunc(2:end-1,1:end-2))/2/dy;
v = -(stmfunc(3:end,2:end-1) - stmfunc(1:end-2,2:end-1))/2/dx;

u = [0, zeros(1,Ny-2), 0;
    zeros(Nx-2,1), u, ones(Nx-2,1) * Unorth(t);
    0, zeros(1,Ny-2), 0];

v = [0, zeros(1,Ny-2), 0;
    zeros(Nx-2,1), v, zeros(Nx-2,1);
    0, zeros(1,Ny-2), 0];

quiver(x2,y2,u',v')
%contourf(x2,y2,stmfunc')
%contourf(x2,y2,vort')



