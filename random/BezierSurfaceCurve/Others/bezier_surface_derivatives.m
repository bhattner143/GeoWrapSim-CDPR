function [B,dBu,dBv,dBu2,dBv2,dBuv] = bezier_surface_derivatives(P,x,y,z)
% Evaluate the Bezier surface and its first and second derivatives at (x,y,z)

% Compute the Bernstein polynomials
u = [1 x x^2 x^3];
v = [1 y y^2 y^3];
Bu = [1 0 0 0; -3 3 0 0; 3 -6 3 0; -1 3 -3 1];
Bv = [1 0 0 0; -3 3 0 0; 3 -6 3 0; -1 3 -3 1];

% Evaluate the surface and its derivatives
P = reshape(P,[4],[4],[3]);

B(:,:,1) = u*P(:,:,1)*Bv'*Bu*v';
B(:,:,2) = u*P(:,:,2)*Bv'*Bu*v';
B(:,:,3) = u*P(:,:,3)*Bv'*Bu*v';

dBu(:,:,1) = [0 1 2*x 3*x^2]*P(:,:,1)*Bv'*Bu*v';
dBu(:,:,2) = [0 1 2*x 3*x^2]*P(:,:,2)*Bv'*Bu*v';
dBu(:,:,3) = [0 1 2*x 3*x^2]*P(:,:,3)*Bv'*Bu*v';

dBv(:,:,1) = u*P(:,:,1)*Bv'*Bu*[0 1 2*y 3*y^2]';
dBv(:,:,2) = u*P(:,:,2)*Bv'*Bu*[0 1 2*y 3*y^2]';
dBv(:,:,3) = u*P(:,:,3)*Bv'*Bu*[0 1 2*y 3*y^2]';

dBu2(:,:,1) = [0 0 2 6*x]*P(:,:,1)*Bv'*Bu*v';
dBu2(:,:,2) = [0 0 2 6*x]*P(:,:,2)*Bv'*Bu*v';
dBu2(:,:,3) = [0 0 2 6*x]*P(:,:,3)*Bv'*Bu*v';

dBv2(:,:,1) = u*P(:,:,1)*Bv'*Bu*[0 0 2 6*y]';
dBv2(:,:,2) = u*P(:,:,2)*Bv'*Bu*[0 0 2 6*y]';
dBv2(:,:,3) = u*P(:,:,3)*Bv'*Bu*[0 0 2 6*y]';

dBuv(:,:,1) = [0 1 2*x 3*x^2]*P(:,:,1)*Bv'*Bu*[0 1 2*y 3*y^2]';
dBuv(:,:,2) = [0 1 2*x 3*x^2]*P(:,:,2)*Bv'*Bu*[0 1 2*y 3*y^2]';
dBuv(:,:,3) = [0 1 2*x 3*x^2]*P(:,:,3)*Bv'*Bu*[0 1 2*y 3*y^2]';