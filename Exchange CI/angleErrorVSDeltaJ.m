dx = 1/2*[0,1;1,0];
dy = 1/2*[0,-1i;1i,0];
dz = 1/2*[1,0;0,-1];
di = 1/2*[1,0;0,1];

kronMat = kron(dx,dx) + kron(dy,dy) + kron(dz,dz);