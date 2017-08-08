function R = quat2rot(q)

if q(1) ~= 1 && q(1)~= -1
    theta = 2 * atan2(sqrt(1 - q(1)^2), q(1));
    k = [q(2) / sqrt(1 - q(1)^2); q(3) / sqrt(1 - q(1)^2); q(4) / sqrt(1 - q(1)^2)];
    v = 1 - cos(theta);
    c = cos(theta);
    s = sin(theta);
    R = [k(1)^2 * v + c,             k(1) * k(2) * v - k(3) * s, k(1) * k(3) * v + k(2) * s;
         k(1) * k(2) * v + k(3) * s, k(2)^2 * v + c,             k(2) * k(3) * v - k(1) * s;
         k(1) * k(3) * v - k(2) * s, k(2) * k(3) * v + k(1) * s, k(3)^2 * v + c];
else
    R = eye(3);
end
end