function RHS = assembleRHS(Nx,Ny,stmfunc,vort,Re,dx,dy,t)
nu = 1/Re;

RHS = zeros(Nx-2,Ny-2);
%interior
for j=2:Ny-3
    for i=2:Nx-3
        fac1 = -(stmfunc(i,j+1) - stmfunc(i,j-1))/2/dy;
        fac2 = (stmfunc(i+1,j) - stmfunc(i-1,j))/2/dx;
        RHS(i,j) = fac1 * (vort(i+1,j)-vort(i-1,j))/2/dx + ...
            fac2 * (vort(i,j+1) - vort(i,j-1))/2/dy + ...
            nu * ((vort(i+1,j) -2*vort(i,j) + vort(i-1,j))/dx/dx+...
                (vort(i,j+1)-2*vort(i,j)+vort(i,j-1))/dy/dy);
    end
end

%South
j = 1;
    for i=2:Nx-3
        vortsouthbc = -2*stmfunc(i,j)/dy/dy;
        fac1 = -(stmfunc(i,j+1) - 0)/2/dy;
        fac2 = (stmfunc(i+1,j) - stmfunc(i-1,j))/2/dx;
        RHS(i,j) = fac1 * (vort(i+1,j)-vort(i-1,j))/2/dx + ...
            fac2 * (vort(i,j+1) - vortsouthbc)/2/dy + ...
            nu * ((vort(i+1,j) -2*vort(i,j) + vort(i-1,j))/dx/dx+...
                (vort(i,j+1)-2*vort(i,j)+vortsouthbc)/dy/dy);
    end
%Southwest
j = 1;
i = 1;
    vortsouthbc = -2*stmfunc(i,j)/dy/dy;
    vortwestbc = -2*stmfunc(i,j)/dx/dx;
    fac1 = -(stmfunc(i,j+1) - 0)/2/dy;
    fac2 = (stmfunc(i+1,j) - 0)/2/dx;
    RHS(i,j) = fac1 * (vort(i+1,j)-vortwestbc)/2/dx + ...
        fac2 * (vort(i,j+1) - vortsouthbc)/2/dy + ...
        nu * ((vort(i+1,j) -2*vort(i,j) + vortwestbc)/dx/dx+...
            (vort(i,j+1)-2*vort(i,j)+vortsouthbc)/dy/dy);
%west
for j=2:Ny-3
    i = 1;
    vortwestbc = -2*stmfunc(i,j)/dx/dx;
    fac1 = -(stmfunc(i,j+1) - stmfunc(i,j-1))/2/dy;
    fac2 = (stmfunc(i+1,j) - 0)/2/dx;
    RHS(i,j) = fac1 * (vort(i+1,j)-vortwestbc)/2/dx + ...
        fac2 * (vort(i,j+1) - vort(i,j-1))/2/dy + ...
        nu * ((vort(i+1,j) -2*vort(i,j) + vortwestbc)/dx/dx+...
            (vort(i,j+1)-2*vort(i,j)+vort(i,j-1))/dy/dy);
end

%Southeast
j = 1;
i = Nx-2;
    vortsouthbc = -2*stmfunc(i,j)/dy/dy;
    vorteastbc = -2*stmfunc(i,j)/dx/dx;
    fac1 = -(stmfunc(i,j+1) - 0)/2/dy;
    fac2 = (0 - stmfunc(i-1,j))/2/dx;
    RHS(i,j) = fac1 * (vorteastbc-vort(i-1,j))/2/dx + ...
        fac2 * (vort(i,j+1) - vortsouthbc)/2/dy + ...
        nu * ((vorteastbc -2*vort(i,j) + vort(i-1,j))/dx/dx+...
            (vort(i,j+1)-2*vort(i,j)+vortsouthbc)/dy/dy);
%East
for j = 2:Ny-3;
    i = Nx-2;
        vorteastbc = -2*stmfunc(i,j)/dx/dx;
        fac1 = -(stmfunc(i,j+1) - stmfunc(i,j-1))/2/dy;
        fac2 = (0 - stmfunc(i-1,j))/2/dx;
        RHS(i,j) = fac1 * (vorteastbc-vort(i-1,j))/2/dx + ...
            fac2 * (vort(i,j+1) - vort(i,j-1))/2/dy + ...
            nu * ((vorteastbc -2*vort(i,j) + vort(i-1,j))/dx/dx+...
                (vort(i,j+1)-2*vort(i,j)+vort(i,j-1))/dy/dy);        
end

%North
j = Ny-2;
    for i = 2:Nx-3
        vortnorthbc = -2*stmfunc(i,j) -Unorth(t)/dy;
        fac1 = -(0 - stmfunc(i,j-1))/2/dy;
        fac2 = (stmfunc(i+1,j) - stmfunc(i-1,j))/2/dx;
        RHS(i,j) = fac1 * (vort(i+1,j)-vort(i-1,j))/2/dx + ...
            fac2 * (vortnorthbc - vort(i,j-1))/2/dy + ...
            nu * ((vort(i+1,j) -2*vort(i,j) + vort(i-1,j))/dx/dx+...
                (vortnorthbc-2*vort(i,j)+vort(i,j-1))/dy/dy);     
    end

%Northwest
j = Ny-2;
    i = 1;
        vortnorthbc = -2*stmfunc(i,j) -Unorth(t)/dy;
        vortwestbc = -2*stmfunc(i,j)/dx/dx;
        fac1 = -(0 - stmfunc(i,j-1))/2/dy;
        fac2 = (stmfunc(i+1,j) - 0)/2/dx;
        RHS(i,j) = fac1 * (vort(i+1,j)-vortwestbc)/2/dx + ...
            fac2 * (vortnorthbc - vort(i,j-1))/2/dy + ...
            nu * ((vort(i+1,j) -2*vort(i,j) + vortwestbc)/dx/dx+...
                (vortnorthbc-2*vort(i,j)+vort(i,j-1))/dy/dy);     
%Northeast
j = Ny-2;
    i = Nx-2;
        vortnorthbc = -2*stmfunc(i,j) -Unorth(t)/dy;
        vorteastbc = -2*stmfunc(i,j)/dx/dx;
        fac1 = -(0 - stmfunc(i,j-1))/2/dy;
        fac2 = (0 - stmfunc(i-1,j))/2/dx;
        RHS(i,j) = fac1 * (vorteastbc-vort(i-1,j))/2/dx + ...
            fac2 * (vortnorthbc - vort(i,j-1))/2/dy + ...
            nu * ((vorteastbc -2*vort(i,j) + vort(i-1,j))/dx/dx+...
                (vortnorthbc-2*vort(i,j)+vort(i,j-1))/dy/dy); 
end