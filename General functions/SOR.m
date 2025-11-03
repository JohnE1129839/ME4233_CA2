function u = SOR(A,b,u0,resobj,omega)

L=tril(A,-1);
D=diag(diag(A));
U=A-L-D;

while 1
    u=(D+omega*L)\(omega*b-(omega*U+(omega-1)*D)*u0);  % SOR method
    res=norm(b-A*u);
    if res<resobj
        break
    end   
    u0=u;
end
end