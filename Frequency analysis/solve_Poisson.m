function [stmfunc] = solve_Poisson(vort,A,Nx,Ny)
    b = -vort(:);
    resobj = 1e-4;
    ngs = (Nx-2) * (Ny-2);
    u0 = zeros(ngs,1);
    stmfunc = SOR(A,b,u0,resobj,1.5);
    stmfunc = reshape(stmfunc, Nx-2, Ny-2);
end