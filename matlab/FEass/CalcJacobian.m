function [detJ, invJt] = Jacobian(xy)
J = [xy(2,1)-xy(1,1) xy(3,1)-xy(1,1); xy(2,2)-xy(1,2) xy(3,2)-xy(1,2)];
detJ = abs(det(J));
invJt = inv(J)';
end