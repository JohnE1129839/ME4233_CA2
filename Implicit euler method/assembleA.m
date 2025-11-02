function A = assembleA(Nx,Ny,dx,dy)
n = (Nx-2) * (Ny-2);
A = zeros(n,n);
c0 = -2/dx^2-2/dy^2;
c1 = 1/dx^2;
c2 = 1/dy^2;
% interior grid points without boundary conditions
for j=2:Ny-3
    for i=2:Nx-2
        po=i+(j-1)*(Nx-2);
        A(po,po)       = c0; % (i,j)
        A(po,po+1)     = c1;   % (i+1,j)
        A(po,po-1)     = c1;   % (i-1,j)
        A(po,po-(Nx-2))= c2;  % (i,j-1)
        A(po,po+(Nx-2))= c2;  % (i,j+1)
    end
end

% South
j=1;
for i=2:Nx-3
        po=i+(j-1)*(Nx-2);
        A(po,po)       = c0; % (i,j)
        A(po,po+1)     = c1;   % (i+1,j)
        A(po,po-1)     = c1;   % (i-1,j)
%        A(po,po-(Nx-1))= c2;  % (i,j-1)
        A(po,po+(Nx-2))= c2;  % (i,j+1)
end

% North
j=Ny-2;
for i=2:Nx-3
        po=i+(j-1)*(Nx-2);
        A(po,po)       = c0; % (i,j)
        A(po,po+1)     = c1;   % (i+1,j)
        A(po,po-1)     = c1;   % (i-1,j)
        A(po,po-(Nx-2))= c2;  % (i,j-1)
%        A(po,po+(Nx-2))= c2;  % (i,j+1)
end

% West
i=1;
for j=2:Ny-3
        po=i+(j-1)*(Nx-2);
        A(po,po)       = c0; % (i,j)
        A(po,po+1)     = c1;   % (i+1,j)
%        A(po,po-1)     = c1;   % (i-1,j)
        A(po,po-(Nx-2))= c2;  % (i,j-1)
        A(po,po+(Nx-2))= c2;  % (i,j+1)
end

% East
i=Nx-2;
for j=2:Ny-3
        po=i+(j-1)*(Nx-2);
        A(po,po)       = c0; % (i,j)
%        A(po,po+1)     = c1;   % (i+1,j)
        A(po,po-1)     = c1;   % (i-1,j)
        A(po,po-(Nx-2))= c2;  % (i,j-1)
        A(po,po+(Nx-2))= c2;  % (i,j+1)
end

% South-west
i=1;j=1;
        po=i+(j-1)*(Nx-2);
        A(po,po)       = c0; % (i,j)
        A(po,po+1)     = c2;   % (i+1,j)
%        A(po,po-1)     = c1;   % (i-1,j)
%        A(po,po-(Nx-2))= c2;  % (i,j-1)
        A(po,po+(Nx-2))= c2;  % (i,j+1)
        

% South-east        
i=Nx-2;j=1;
        po=i+(j-1)*(Nx-2);
        A(po,po)       = c0; % (i,j)
%        A(po,po+1)     = c1;   % (i+1,j)
        A(po,po-1)     = c1;   % (i-1,j)
%        A(po,po-(Nx-2))= c2;  % (i,j-1)
        A(po,po+(Nx-2))= c2;  % (i,j+1)

% North-east         
i=Nx-2;j=Ny-2;
        po=i+(j-1)*(Nx-2);
        A(po,po)       = c0; % (i,j)
%        A(po,po+1)     = c1;   % (i+1,j)
        A(po,po-1)     = c1;   % (i-1,j)
        A(po,po-(Nx-2))= c2;  % (i,j-1)
%        A(po,po+(Nx-2))= c2;  % (i,j+1)

% North-west           
i=1;j=Ny-2;
        po=i+(j-1)*(Nx-2);
        A(po,po)       = c0; % (i,j)
        A(po,po+1)     = c1;   % (i+1,j)
%        A(po,po-1)     = c1;   % (i-1,j)
        A(po,po-(Nx-2))= c2;  % (i,j-1)
%        A(po,po+(Nx-2))= c2;  % (i,j+1)
        
end