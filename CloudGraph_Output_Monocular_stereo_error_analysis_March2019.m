clear all;
%%
%Athroscope positions (current = A2, previous = A1).
%Assumes that the rotation has been accounted for (i.e. pure translation).
%0 Deg
% A1 = [0 0 0]
% A2 = [0 0 1]

%45 Deg
A1 = [0 0 0]
A2 = [-2 2 -0.7]
% A2 = [-4 4 -1.4]

%90 Deg
% A1 = [0 0 0]
% A2 = [3 0 -2]

% Select two points to compute the distance (e.g. two points detected on the knee gap).
G1 = [-1 3 0]
G2 = [1 2 -1]

%plot positions
figure(1)
scatter3(A1(1), A1(2), A1(3), 'b')
hold on;
scatter3(A2(1), A2(2), A2(3), 'b')
scatter3(G1(1), G1(2), G1(3), 'r')
scatter3(G2(1), G2(2), G2(3), 'r')
A = [A1;A2]
G = [G1;G2]
plot3(A(:,1),A(:,2),A(:,3),'green')
plot3(G(:,1),G(:,2),G(:,3),'yellow')

%Compute the ground truth distance
d_gt = pdist([G1;G2]);
%Translation of athroscope tip
T = A2 - A1

%Compute the vectors which will be required for the angle computations.
%The actual vectors will be unknon in the knee gap but the unit vectors
%can be determined using the azimuth and elevation computed by a camera
%model. 
%Azimuth and elevation to unit vector:
% https://math.stackexchange.com/questions/1150232/finding-the-unit-direction-vector-given-azimuth-and-elevation
a1 = G1 - A2
a2 = G2 - A2
b1 = G1 - A1;
b2 = G2 - A1;

%Compute the angles we know
alpha1 = acos(dot(T/norm(T),b1/norm(b1)));
alpha2 = acos(dot(T/norm(T),b2/norm(b2)));
beta1 = acos(-dot(T/norm(T),a1/norm(a1)));
beta2 = acos(-dot(T/norm(T),a2/norm(a2)));

%gamma will be useful for the sine rule in the next step
gamma1 = pi()-alpha1-beta1;
gamma2 = pi()-alpha2-beta2;

%using the sine rule to compute the length of a1 and a2
a1_len = sin(alpha1)/sin(gamma1) * norm(T);
a2_len = sin(alpha2)/sin(gamma2) * norm(T);
b1_len = sin(beta1)/sin(gamma1) * norm(T);
b2_len = sin(beta2)/sin(gamma2) * norm(T);

%now compute the a1 and a2 vectors
a1_computed = a1_len * (a1/norm(a1))
a2_computed = a2_len * (a2/norm(a2))
b1_computed = b1_len * (b1/norm(b1))
b2_computed = b2_len * (b2/norm(b2))

%Compute the distance between a1 and a2
d_gt
d = pdist([a2_computed;a1_computed])
error = d - d_gt

%% Initial error setup
figure(1)

dT = 0.0367; %RMS Value of translation error
image_error = deg2rad(2.367);%Set 3 Image Error
% image_error = deg2rad(4.302);%Set 2 Image Error
% image_error = deg2rad(5.586);%Set 6 Image Error

o = [0 0 0]
scatter3(o(1),o(2),o(3))
qiver_T = quiver3(o(1),o(2),o(3),T(1),T(2),T(3),0)
qiver_T.Color = 'green'
qiver_a1 = quiver3(T(1),T(2),T(3),a1(1),a1(2),a1(3),0)
qiver_a1.Color = 'blue'
qiver_a2 = quiver3(T(1),T(2),T(3),a2(1),a2(2),a2(3),0)
qiver_a2.Color = 'cyan'
qiver_b1 = quiver3(o(1),o(2),o(3),b1(1),b1(2),b1(3),0)
qiver_b1.Color = [0.5,0.5,0.5]
qiver_b2 = quiver3(o(1),o(2),o(3),b2(1),b2(2),b2(3),0)
qiver_b2.Color = 'red'

%% Translation error volumes
%Define the number of points in the error volume/cloud
numStepsEl = 4
numStepsAz = 8
azStepSize = deg2rad(360)/numStepsAz
elStepSize = deg2rad(180)/numStepsEl
[az_1,el_1]=meshgrid((deg2rad(0):azStepSize:deg2rad(360)),(deg2rad(-90):elStepSize:deg2rad(90)));

%Construct translation error volumes (spheres)
Rx_1 = sin(az_1).*cos(el_1);
Ry_1 = cos(az_1).*cos(el_1);
Rz_1 = sin(el_1);
Sx_1 = Rx_1*dT;
Sy_1 = Ry_1*dT;
Sz_1 = Rz_1*dT;

[az_2,el_2]=meshgrid((deg2rad(0):azStepSize:deg2rad(360)),(deg2rad(-90):elStepSize:deg2rad(90)));
Rx_2 = sin(az_2).*cos(el_2);
Ry_2 = cos(az_2).*cos(el_2);
Rz_2 = sin(el_2);
Sx_2 = T(1) + Rx_2*dT;
Sy_2 = T(2) + Ry_2*dT;
Sz_2 = T(3) + Rz_2*dT;

%Reformat data for vector based computation (much faster!)
SX_1 = Sx_1(:);
SY_1 = Sy_1(:);
SZ_1 = Sz_1(:);

SX_2 = Sx_2(:);
SY_2 = Sy_2(:);
SZ_2 = Sz_2(:);

%Compute the various potential translation directions due to error
[i1,i2] = meshgrid((1:length(SX_1)),(1:length(SX_2)));
T_dash_X = SX_2(i2)-SX_1(i1);
T_dash_Y = SY_2(i2)-SY_1(i1);
T_dash_Z = SZ_2(i2)-SZ_1(i1);
T_dash_X = T_dash_X(:);
T_dash_Y = T_dash_Y(:);
T_dash_Z = T_dash_Z(:);

%Match translation vectors to error volume points
SX_1_array = SX_1(i1);
SY_1_array = SY_1(i1);
SZ_1_array = SZ_1(i1);

SX_1_array = SX_1_array(:);
SY_1_array = SY_1_array(:);
SZ_1_array = SZ_1_array(:);

SX_2_array = SX_2(i1);
SY_2_array = SY_2(i1);
SZ_2_array = SZ_2(i1);

SX_2_array = SX_2_array(:);
SY_2_array = SY_2_array(:);
SZ_2_array = SZ_2_array(:);

%% Create distance error volume
numSteps = 4
coneStepSize = deg2rad(360)/numSteps
T_vec = [T(1); T(2); T(3)];

for i=1:length(T_dash_X)
    %Reformate the tranlation vectors
    T_dash = [T_dash_X(i); T_dash_Y(i);T_dash_Z(i)];
    S1 = [SX_1_array(i); SY_1_array(i);SZ_1_array(i)];
    S2 = [SX_2_array(i); SY_2_array(i);SZ_2_array(i)];
    T_dash_norm = T_dash/norm(T_dash);
    %Compute dAlpha, dBeta, dGama, dL
    %Compute the new angles due to errors
    alpha1_dash = acos(dot(T_dash_norm,b1/norm(b1)));
    alpha2_dash = acos(dot(T_dash_norm,b2/norm(b2)));
    beta1_dash = acos(-dot(T_dash_norm,a1/norm(a1)));
    beta2_dash = acos(-dot(T_dash_norm,a2/norm(a2)));

    %gamma will be useful for the sine rule in the next step
    gamma1_dash = pi()-alpha1_dash-beta1_dash;
    gamma2_dash = pi()-alpha2_dash-beta2_dash;

    %using the sine rule to compute the length of a1 and a2
    a1_len_dash = sin(alpha1_dash)/sin(gamma1_dash) * norm(T_dash);
    a2_len_dash = sin(alpha2_dash)/sin(gamma2_dash) * norm(T_dash);    
    b1_len_dash = sin(beta1_dash)/sin(gamma1_dash) * norm(T_dash);
    b2_len_dash = sin(beta2_dash)/sin(gamma2_dash) * norm(T_dash); 
    
    %Add image cone error b1
    %a1
    a1_norm = a1/norm(a1);
    dL_a1 =  a1_len - a1_len_dash;
    L_a1 = a1_len;
    dG_a1 = (cross(a1_norm,cross(T_dash_norm,a1_norm)))/norm(cross(a1_norm,cross(T_dash_norm,a1_norm)))*sin(image_error)*(L_a1+dL_a1);
    dG_img_error_a1 = a1_norm*(L_a1+dL_a1) + dG_a1;

    %a2
    a2_norm = a2/norm(a2);
    dL_a2 =  a2_len - a2_len_dash;
    L_a2 = a2_len;
    dG_a2 = (cross(a2_norm,cross(T_dash_norm,a2_norm)))/norm(cross(a2_norm,cross(T_dash_norm,a2_norm)))*sin(image_error)*(L_a2+dL_a2);
    dG_img_error_a2 = a2_norm*(L_a2+dL_a2) + dG_a2;
    %b1
    b1_norm = b1/norm(b1);
    dL_b1 =  b1_len - b1_len_dash;
    L_b1 = b1_len;
    dG_b1 = (cross(b1_norm,cross(T_dash_norm,b1_norm)))/norm(cross(b1_norm,cross(T_dash_norm,b1_norm)))*sin(image_error)*(L_b1+dL_b1);
    dG_img_error_b1 = b1_norm*(L_b1+dL_b1) + dG_b1;
    %b2
    b2_norm = b2/norm(b2);
    dL_b2 =  b2_len - b2_len_dash;
    L_b2 = b2_len;
    dG_b2 = (cross(b2_norm,cross(T_dash_norm,b2_norm)))/norm(cross(b2_norm,cross(T_dash_norm,b2_norm)))*sin(image_error)*(L_b2+dL_b2);
    dG_img_error_b2 = b2_norm*(L_b2+dL_b2) + dG_b2;
    j = 1;
    for theta = deg2rad(0):coneStepSize:deg2rad(360)
        %a1
        v_rot_a1 = transpose(S2)+rodrigues_rot(dG_img_error_a1,a1_norm,theta);
        v_rot_x_a1(i,j) = v_rot_a1(1);
        v_rot_y_a1(i,j) = v_rot_a1(2);
        v_rot_z_a1(i,j) = v_rot_a1(3);
        %a2
        v_rot_a2 = transpose(S2)+rodrigues_rot(dG_img_error_a2,a2_norm,theta);
        v_rot_x_a2(i,j) = v_rot_a2(1);
        v_rot_y_a2(i,j) = v_rot_a2(2);
        v_rot_z_a2(i,j) = v_rot_a2(3);
        %b1
        v_rot_b1 = transpose(S1)+rodrigues_rot(dG_img_error_b1,b1_norm,theta);
        v_rot_x_b1(i,j) = v_rot_b1(1);
        v_rot_y_b1(i,j) = v_rot_b1(2);
        v_rot_z_b1(i,j) = v_rot_b1(3);
        %b2
        v_rot_b2 = transpose(S1)+rodrigues_rot(dG_img_error_b2,b2_norm,theta);
        v_rot_x_b2(i,j) = v_rot_b2(1);
        v_rot_y_b2(i,j) = v_rot_b2(2);
        v_rot_z_b2(i,j) = v_rot_b2(3);
        j = j + 1;
    end

end
%a1
error_vol_x_a1 = v_rot_x_a1(:);
error_vol_y_a1 = v_rot_y_a1(:);
error_vol_z_a1 = v_rot_z_a1(:);
%a2
error_vol_x_a2 = v_rot_x_a2(:);
error_vol_y_a2 = v_rot_y_a2(:);
error_vol_z_a2 = v_rot_z_a2(:);
%b1
error_vol_x_b1 = v_rot_x_b1(:);
error_vol_y_b1 = v_rot_y_b1(:);
error_vol_z_b1 = v_rot_z_b1(:);
%b2
error_vol_x_b2 = v_rot_x_b2(:);
error_vol_y_b2 = v_rot_y_b2(:);
error_vol_z_b2 = v_rot_z_b2(:);
%% Plotting
figure(1)
mesh(Sx_1,Sy_1,Sz_1)
hold on;
mesh(Sx_2,Sy_2,Sz_2)
%a1
scatter3(error_vol_x_a1,error_vol_y_a1,error_vol_z_a1, 'MarkerEdgeColor','b','MarkerFaceColor', 'b')
%a2
scatter3(error_vol_x_a2,error_vol_y_a2,error_vol_z_a2, 'MarkerEdgeColor','cyan','MarkerFaceColor', 'cyan')
%b1
scatter3(error_vol_x_b1,error_vol_y_b1,error_vol_z_b1, 'MarkerEdgeColor',[0.5,0.5,0.5],'MarkerFaceColor', [0.5,0.5,0.5])
%b2
scatter3(error_vol_x_b2,error_vol_y_b2,error_vol_z_b2, 'MarkerEdgeColor','black','MarkerFaceColor', 'red')
% DT = delaunayTriangulation(x,y,z);  %# Create the tetrahedral mesh
% tetramesh(DT);

axis equal
% xyzlabel()
hold off;

%%
%%Compute the max and minimum distances from the point clouds.
% error_a1 = cat(2,error_vol_x_a1,error_vol_y_a1,error_vol_z_a1);
% error_a2 = cat(2,error_vol_x_a2,error_vol_y_a2,error_vol_z_a2);
% error_b1 = cat(2,error_vol_x_b1,error_vol_y_b1,error_vol_z_b1);
% error_b2 = cat(2,error_vol_x_b2,error_vol_y_b2,error_vol_z_b2);
% 
% error_ab1 = [error_a1;error_b1];
% size(error_ab1)
% error_ab2 = [error_a2;error_b2];
% size(error_ab2)
% 
% [ii1,ii2] = meshgrid((1:length(error_ab1)),(1:length(error_ab2)));
% d_cloud_array = pdist([error_ab1(ii1);error_ab2(ii2)]);
% size(d_cloud_array)
% d_max = max(d_cloud_array);
% d_min = min(d_cloud_array);

%% Compute the sphere approximation
tic
theta_T = atan(2*dT/norm(T));

ra1_T = (L_a1+dL_a1)*sin(theta_T);
ra2_T = (L_a2+dL_a2)*sin(theta_T);
rb1_T = (L_b1+dL_b1)*sin(theta_T);
rb2_T = (L_b2+dL_b2)*sin(theta_T);

ra1_psi = (L_a1+dL_a1)*sin(image_error);
ra2_psi = (L_a2+dL_a2)*sin(image_error);
rb1_psi = (L_b1+dL_b1)*sin(image_error);
rb2_psi = (L_b2+dL_b2)*sin(image_error);

R_a1 = ra1_T + ra1_psi;
R_a2 = ra2_T + ra2_psi;
R_b1 = rb1_T + rb1_psi;
R_b2 = rb2_T + rb2_psi;

[x_a1,y_a1,z_a1] = sphere_points(R_a1,T(1) + a1_computed(1),T(2) + a1_computed(2),T(3) + a1_computed(3));
[x_a2,y_a2,z_a2] = sphere_points(R_a2,T(1) + a2_computed(1),T(2) + a2_computed(2),T(3) + a2_computed(3));
[x_b1,y_b1,z_b1] = sphere_points(R_b1,b1_computed(1),b1_computed(2),b1_computed(3));
[x_b2,y_b2,z_b2] = sphere_points(R_b2,b2_computed(1),b2_computed(2),b2_computed(3));


d_a_min = norm(a2_computed-a1_computed)-R_a1-R_a2;
d_a_max = norm(a2_computed-a1_computed)+R_a1+R_a2;
d_b_min = norm(b2_computed-b1_computed)-R_b1-R_b2;
d_b_max = norm(b2_computed-b1_computed)+R_b1+R_b2;

d_gt
d_min_sphere_approx = min([d_a_min,d_b_min]);
d_max_sphere_approx = max([d_a_max,d_b_max]);
range_sphere_approx = [d_min_sphere_approx, d_max_sphere_approx]
toc
%Plot spheres
figure(1);
hold on;
s = surf(x_a1,y_a1,z_a1,'FaceAlpha',0.2);
s.EdgeColor = [0.2,0.2,0.2];
s.FaceColor = 'blue';
s = surf(x_a2,y_a2,z_a2,'FaceAlpha',0.2);
s.EdgeColor = [0.2,0.2,0.2];
s.FaceColor = 'cyan';
s = surf(x_b1,y_b1,z_b1,'FaceAlpha',0.2);
s.EdgeColor = [0.2,0.2,0.2];
s.FaceColor = [0.5,0.5,0.5];
s = surf(x_b2,y_b2,z_b2,'FaceAlpha',0.2);
s.EdgeColor = [0.2,0.2,0.2];
s.FaceColor = 'red';




