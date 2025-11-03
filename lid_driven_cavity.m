timeUpdateMethod = "I"; % E = Explicit, I = implicit

%Importing general functions
addpath("./General functions");

%Initializing parameters
T = readtable("./Parameters.txt");
Re = T.Re;
Nx = T.Nx; Ny = T.Ny;
Lx = T.Lx; Ly = T.Ly;

x = linspace(0,Lx,Nx); dx = x(2) -x(1);
y = linspace(0,Ly,Ny); dy = y(2) - y(1);

[x2,y2] = meshgrid(x,y);
if timeUpdateMethod == "E"
    dt = T.dte; 
else
    dt = T.dti;
end

tf = T.tf; 
Nmax = tf/dt;

%Finding indices of point c
ic = T.cx/dx; jc = T.cy/dy;
c_u = [];

%Generating A matrix
[A] = assembleA(Nx,Ny,dx,dy);

%Initializing vort matrix
vort = zeros(Nx-2,Ny-2);
t = 0;

%Solving for streamfunction
for i = 1:Nmax
    stmfunc = solve_Poisson(vort,A,Nx,Ny);
    if timeUpdateMethod == "E"
        vort = advanceVortExplicit(stmfunc,vort,Nx,Ny,dx,dy,dt,Re,t);
    else
        vort = advanceVortImplicit(stmfunc,vort,Nx,Ny,dx,dy,dt,Re,t);
    end
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

figure
quiver(x2,y2,u',v')
title("Flow velocity at tf");
figure
contourf(x2,y2,stmfunc')
title("Streamfunction at tf");

%% 

%plotting u on c
figure;
plot((1:Nmax)*dt, c_u);
xlabel('Time (s)');
ylabel('Amplitude');
title('Time Domain Signal: c\_u(t)');
grid on;
figure;
%Fourrier analysis on c
fs = 1/dt;
N = length(c_u);
Y = fft(c_u);
P2 = abs(Y/N);
P1 = P2(1:N/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = fs * (0:(N/2)) / N;
plot(f, P1);
xlabel('Frequency (Hz)');
ylabel('Amplitude');
title('Single-Sided Amplitude Spectrum of c\_u(t)');
grid on;

