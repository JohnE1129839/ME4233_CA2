function vortnew = advanceVortImplicit(stmfunc,vort,Nx,Ny,dx,dy,dt,Re,t)
%Using the backward euler discretization
C = dx^2/dt; k = dx/dy;
ngs = (Nx-2) * (Ny-2);
A = zeros(ngs);
vort = reshape(vort,ngs,1);
b(:) = 4 * C * vort;

%Interior points
for i=2:Nx-3
    for j=2:Ny-3
        po = i + (j-1) * (Nx-2);
        A(po,po)       = (4 * C + 8/Re *(1+k^2)); % (i,j)
        A(po,po+1)     = k*(stmfunc(i,j+1)-stmfunc(i,j-1)) - 4/Re;% (i+1,j)
        A(po,po-1)     = -k*(stmfunc(i,j+1)-stmfunc(i,j-1)) - 4/Re;% (i-1,j)
        A(po,po+(Nx-2))= -k*(stmfunc(i+1,j)-stmfunc(i-1,j)) - k^2 * 4/Re;  % (i,j-1)
        A(po,po-(Nx-2))= k*(stmfunc(i+1,j)-stmfunc(i-1,j)) - k^2 * 4/Re;  % (i,j+1)
    end
end

%North points
for i=2:Nx-3
    j=Ny-2;
        vortnorthbc = -2*stmfunc(i,j)/dy/dy -Unorth(t,(i-1)*dx,(j-1)*dy)/dy;
        po = i + (j-1) * (Nx-2);
        A(po,po)       = (4 * C + 8/Re *(1+k^2)); % (i,j)
        A(po,po+1)     = k*(0-stmfunc(i,j-1)) - 4/Re;% (i+1,j)
        A(po,po-1)     = -k*(0-stmfunc(i,j-1)) - 4/Re;% (i-1,j)
        %A(po,po+(Nx-2))= -k(stmfunc(i+1,j)-stmfunc(i-1,j)) - k^2 * 4/Re;  % (i,j+1)
        A(po,po-(Nx-2))= k*(stmfunc(i+1,j)-stmfunc(i-1,j)) - k^2 * 4/Re;  % (i,j-1)
    b(po) = b(po) - (-k*(stmfunc(i+1,j)-stmfunc(i-1,j)) - k^2 * 4/Re) * (vortnorthbc);
end

%South points
for i=2:Nx-3
    j = 1;
        vortsouthbc = -2*stmfunc(i,j)/dy/dy;
        po = i + (j-1) * (Nx-2);
        A(po,po)       = (4 * C + 8/Re *(1+k^2)); % (i,j)
        A(po,po+1)     = k*(stmfunc(i,j+1)-0) - 4/Re;% (i+1,j)
        A(po,po-1)     = -k*(stmfunc(i,j+1)-0) - 4/Re;% (i-1,j)
        A(po,po+(Nx-2))= -k*(stmfunc(i+1,j)-stmfunc(i-1,j)) - k^2 * 4/Re;  % (i,j+1)
        %A(po,po-(Nx-2))= k(stmfunc(i+1,j)-stmfunc(i-1,j)) - k^2 * 4/Re;  % (i,j-1)
        b(po) = b(po) - (k*(stmfunc(i+1,j)-stmfunc(i-1,j)) - k^2 * 4/Re) * vortsouthbc;
end

%East points
i = Nx-2;
    for j=2:Ny-3
        vorteastbc = -2*stmfunc(i,j)/dx/dx;
        po = i + (j-1) * (Nx-2);
        A(po,po)       = (4 * C + 8/Re *(1+k^2)); % (i,j)
        %A(po,po+1)     = k(stmfunc(i,j+1)-stmfunc(i,j-1)) - 4/Re;% (i+1,j)
        A(po,po-1)     = -k*(stmfunc(i,j+1)-stmfunc(i,j-1)) - 4/Re;% (i-1,j)
        A(po,po+(Nx-2))= -k*(0-stmfunc(i-1,j)) - k^2 * 4/Re;  % (i,j+1)
        A(po,po-(Nx-2))= k*(0-stmfunc(i-1,j)) - k^2 * 4/Re;  % (i,j-1)
        b(po) = b(po) - (k*(stmfunc(i,j+1)-stmfunc(i,j-1)) - 4/Re) * vorteastbc;
    end

%West points
i = 1;
    for j=2:Ny-3
        vortwestbc = -2*stmfunc(i,j)/dx/dx;
        po = i + (j-1) * (Nx-2);
        A(po,po)       = (4 * C + 8/Re *(1+k^2)); % (i,j)
        A(po,po+1)     = k*(stmfunc(i,j+1)-stmfunc(i,j-1)) - 4/Re;% (i+1,j)
        %A(po,po-1)     = -k*(stmfunc(i,j+1)-stmfunc(i,j-1)) - 4/Re;% (i-1,j)
        A(po,po+(Nx-2))= -k*(stmfunc(i+1,j)-0) - k^2 * 4/Re;  % (i,j+1)
        A(po,po-(Nx-2))= k*(stmfunc(i+1,j)-0) - k^2 * 4/Re;  % (i,j-1)
        b(po) = b(po) - (-k*(stmfunc(i,j+1)-stmfunc(i,j-1)) - 4/Re) * vortwestbc;
    end

%NorthEast points
i = Nx-2;
    j=Ny-2;
        vortnorthbc = -2*stmfunc(i,j)/dy/dy -Unorth(t,(i-1)*dx,(j-1)*dy)/dy;
        vorteastbc = -2*stmfunc(i,j)/dx/dx;
        po = i + (j-1) * (Nx-2);
        A(po,po)       = (4 * C + 8/Re *(1+k^2)); % (i,j)
        %A(po,po+1)     = k(stmfunc(i,j+1)-stmfunc(i,j-1)) - 4/Re;% (i+1,j)
        A(po,po-1)     = -k*(0-stmfunc(i,j-1)) - 4/Re;% (i-1,j)
        %A(po,po+(Nx-2))= -k(stmfunc(i+1,j)-stmfunc(i-1,j)) - k^2 * 4/Re;  % (i,j+1)
        A(po,po-(Nx-2))= k*(0-stmfunc(i-1,j)) - k^2 * 4/Re;  % (i,j-1)
    b(po) = b(po) - (-k*(0-stmfunc(i-1,j)) - k^2 * 4/Re) * (vortnorthbc)...
        -(k*(0-stmfunc(i,j-1)) - 4/Re) * vorteastbc;

%NorthWest points
i = 1;
    j=Ny-2;
        vortnorthbc = -2*stmfunc(i,j)/dy/dy -Unorth(t,(i-1)*dx,(j-1)*dy)/dy;
        vortwestbc = -2*stmfunc(i,j)/dx/dx;
        po = i + (j-1) * (Nx-2);
        A(po,po)       = (4 * C + 8/Re *(1+k^2)); % (i,j)
        A(po,po+1)     = k*(0-stmfunc(i,j-1)) - 4/Re;% (i+1,j)
        %A(po,po-1)     = -k(stmfunc(i,j+1)-stmfunc(i,j-1)) - 4/Re;% (i-1,j)
        %A(po,po+(Nx-2))= -k(stmfunc(i+1,j)-stmfunc(i-1,j)) - k^2 * 4/Re;  % (i,j+1)
        A(po,po-(Nx-2))= k*(stmfunc(i+1,j)-0) - k^2 * 4/Re;  % (i,j-1)
    b(po) = b(po) - (-k*(stmfunc(i+1,j)-0) - k^2 * 4/Re) * (vortnorthbc)...
        -(-k*(0-stmfunc(i,j-1)) - 4/Re) * vortwestbc;

%SouthWest points
i = 1;
    j = 1;
        vortwestbc = -2*stmfunc(i,j)/dx/dx;
        vortsouthbc = -2*stmfunc(i,j)/dy/dy;
        po = i + (j-1) * (Nx-2);
        A(po,po)       = (4 * C + 8/Re *(1+k^2)); % (i,j)
        A(po,po+1)     = k*(stmfunc(i,j+1)-0) - 4/Re;% (i+1,j)
        %A(po,po-1)     = -k(stmfunc(i,j+1)-stmfunc(i,j-1)) - 4/Re;% (i-1,j)
        A(po,po+(Nx-2))= -k*(stmfunc(i+1,j)-0) - k^2 * 4/Re;  % (i,j+1)
        %A(po,po-(Nx-2))= k(stmfunc(i+1,j)-stmfunc(i-1,j)) - k^2 * 4/Re;  % (i,j-1)
        b(po) = b(po) - (k*(stmfunc(i+1,j)-0) - k^2 * 4/Re) * vortsouthbc ...
            -(-k*(stmfunc(i,j+1)-0) - 4/Re) * vortwestbc;

%SouthEastpoints
i = Nx-2;
    j = 1;
        vorteastbc = -2*stmfunc(i,j)/dx/dx;
        vortsouthbc = -2*stmfunc(i,j)/dy/dy;
        po = i + (j-1) * (Nx-2);
        A(po,po)       = (4 * C + 8/Re *(1+k^2)); % (i,j)
        %A(po,po+1)     = k(stmfunc(i,j+1)-stmfunc(i,j-1)) - 4/Re;% (i+1,j)
        A(po,po-1)     = -k*(stmfunc(i,j+1)-0) - 4/Re;% (i-1,j)
        A(po,po+(Nx-2))= -k*(0-stmfunc(i-1,j)) - k^2 * 4/Re;  % (i,j+1)
        %A(po,po-(Nx-2))= k(stmfunc(i+1,j)-stmfunc(i-1,j)) - k^2 * 4/Re;  % (i,j-1)
        b(po) = b(po) - (k*(0-stmfunc(i-1,j)) - k^2 * 4/Re) * vortsouthbc ...
            -(k*(stmfunc(i,j+1)-0) - 4/Re) * vorteastbc;
vortnew = A\b';       
vortnew = reshape(vortnew,Nx-2,Ny-2);
end