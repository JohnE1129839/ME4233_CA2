function vortnew = advanceVortExplicit(stmfunc,vort,Nx,Ny,dx,dy,dt,Re,t)
%Using the forward euler time discretization
RHS = assembleRHS(Nx,Ny,stmfunc,vort,Re,dx,dy,t);
vortnew = vort + dt * RHS;
end