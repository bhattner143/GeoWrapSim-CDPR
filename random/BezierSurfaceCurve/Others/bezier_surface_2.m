function B = bezier_surface_2(P, n, u, v, w)
% Evaluate the Bezier surface defined by control points P at (u,v,w)

% Compute the Bernstein polynomials for u, v, and w
Bu = bernstein(n, u);
Bv = bernstein(n, v);
Bw = bernstein(n, w);

P_array(:,:,1) = reshape(P(:,1),4,[]);
P_array(:,:,2) = reshape(P(:,2),4,[]);
P_array(:,:,3) = reshape(P(:,3),4,[]);

% Compute the weighted sum of the control points
B = zeros(1,1,1);
for i = 1:3
    B(:,:,i) = Bu'*P_array(:,:,1)*Bv;
end
B = reshape(B,[3],[1]);
% for i = 0:n
%     for j = 0:n
%         for k = 0:n
%             B = B + (P(i*(n+1)^2 + j*(n+1) + k+1,:) * ...
%                 Bu(i+1) * Bv(j+1) * Bw(k+1)).';
%         end
%     end
% end
end