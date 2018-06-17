function dz=reducedbifunc(z,G,M,S,L,v)
w=z(1:v);
x_r=z(v+1:end);

dw=S*w;
dx_r=(S-G*L-M*x_r*L)*x_r+G*L*w+M*x_r*L*w;
dz=[dw;dx_r];
end
