%Initializing parameters
Re = 27;
Nx = 51; Ny = 21;
Lx = 3; Ly = 1;

x = linspace(0,Lx,Nx); dx = x(2) -x(1);
y = linspace(0,Ly,Ny); dy = y(2) - y(1);
ic = 1.5/dx; jc = 0.5/dy;

[x2,y2] = meshgrid(x,y);
c_u = [];

dt = 0.01; tf = 1; Nmax = tf/dt;

%Generating A matrix
[A] = assembleA(Nx,Ny,dx,dy);

%Initializing vort matrix
vort = zeros(Nx-2,Ny-2);
t = 0;

%Solving for streamfunction
for i = 1:Nmax
    stmfunc = solve_Poisson(vort,A,Nx,Ny);
    vort = advanceVort(stmfunc,vort,Nx,Ny,dx,dy,dt,Re,t);
    c_u = [c_u (stmfunc(ic,jc+1)-stmfunc(ic,jc-1))/2/dy];
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



