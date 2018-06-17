function dz=bifunc(z,A,B,N,S,L,v)
w=z(1:v);
x=z(v+1:end);

dw=S*w;
dx=A*x+B*L*w+N*x*L*w;
dz=[dw;dx];
end
