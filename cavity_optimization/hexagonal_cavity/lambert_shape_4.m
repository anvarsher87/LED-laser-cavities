% lambertian
function dir_s=lambert_shape_4 (n, tang1, tang2, solid_angle)
sin_theta = (rand).^(1/2).*sin(solid_angle);

cos_theta = sqrt(1-sin_theta.*sin_theta);

%random in plane angle
psi = rand*2*pi;

% three vector components
% a = sin_theta.*cos(psi);
% b = sin_theta.*sin(psi);
% c = cos_theta;

a = sin_theta.*cos(psi);
b = sin_theta.*sin(psi);
c = cos_theta;

% replicate basis vectors along second dimension
tang1_rep = repmat(tang1,  1);
tang2_rep = repmat(tang2,  1);
n_rep = repmat(n,  1);

% calculate vector components
v1 = a' .* tang1_rep;
v2 = b' .* tang2_rep;
v3 = c' .* n_rep;

dir_s = v1+v2+v3;
