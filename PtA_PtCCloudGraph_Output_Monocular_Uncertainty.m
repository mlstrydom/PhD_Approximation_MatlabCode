clear all;
close all;
compute_sphere_approx = 1;
compute_error_cloud = 1;
compute_cloud_range = 1*compute_error_cloud;
plot_cloud_points = 1;
%%
%Translation errors
numStepsEl = 4;
numStepsAz = 6;
%image errors
numSteps = 8 ;
% Translation angle
T_Angle=0;
% Translation Distance
T_dist = 2;
% Translation Increments
dT = 0.0367 ;%RMS Value of translation error
%dT = 0.0698 %RMS Value of translation error
%Set Gap size [mm]
Gap_size = 4;

%%
%Athroscope positions (current = A2, previous = A1).
%Assumes that the rotation has been accounted for (i.e. pure translation).
dir = deg2rad(T_Angle);
%Athroscope positions (current = A2, previous = A1).
%Assumes that the rotation has been accounted for (i.e. pure translation).
A1_tmp = [-T_dist*sin(dir) -T_dist*cos(dir) 0];
A1 = [-T_dist*sin(dir) -T_dist*cos(dir) 0]-A1_tmp;
A2 = [0 0 0]-A1_tmp;

%Select two points to compute the distance (e.g. two points detected on the knee gap).
G1 = [-1 3 0];
G2 = [1 2 1];
G_dir = (G2-G1)/norm(G2-G1);
G2 = G1+Gap_size*G_dir; %Scale the gap in the direction

G1 = G1-A1_tmp;
G2 = G2-A1_tmp;

%plot positions
figure(1)
scatter3(A1(1), A1(2), A1(3), 'b');
hold on;
scatter3(A2(1), A2(2), A2(3), 'b');
scatter3(G1(1), G1(2), G1(3), 'r');
scatter3(G2(1), G2(2), G2(3), 'r');
A = [A1;A2];
G = [G1;G2];
plot3(A(:,1),A(:,2),A(:,3),'green');
plot3(G(:,1),G(:,2),G(:,3),'yellow');
view([-15 -30 20])

%Compute the ground truth distance
d_gt = pdist([G1;G2]);
%Translation of athroscope tip
T = A2 - A1;
T_norm = T/norm(T);


%Compute the vectors which will be required for the angle computations.
%The actual vectors will be unknon in the knee gap but the unit vectors
%can be determined using the azimuth and elevation computed by a camera
%model. 
%Azimuth and elevation to unit vector:
% https://math.stackexchange.com/questions/1150232/finding-the-unit-direction-vector-given-azimuth-and-elevation
a1 = G1 - A2;
a2 = G2 - A2;
b1 = G1 - A1;
b2 = G2 - A1;

%Compute the angles we know
alpha1 = acos(dot(T_norm,b1/norm(b1)));
alpha2 = acos(dot(T_norm,b2/norm(b2)));
beta1 = acos(-dot(T_norm,a1/norm(a1)));
beta2 = acos(-dot(T_norm,a2/norm(a2)));

%gamma will be useful for the sine rule in the next step
gamma1 = pi()-alpha1-beta1;
gamma2 = pi()-alpha2-beta2;

%using the sine rule to compute the length of a1 and a2
a1_len = sin(alpha1)/sin(gamma1) * norm(T);
a2_len = sin(alpha2)/sin(gamma2) * norm(T);
b1_len = sin(beta1)/sin(gamma1) * norm(T);
b2_len = sin(beta2)/sin(gamma2) * norm(T);

L_a1 = a1_len;
L_a2 = a2_len;
L_b1 = b1_len;
L_b2 = b2_len;

%now compute the a1 and a2 vectors
a1_computed = a1_len * (a1/norm(a1));
a2_computed = a2_len * (a2/norm(a2));
b1_computed = b1_len * (b1/norm(b1));
b2_computed = b2_len * (b2/norm(b2));

%Compute the distance between a1 and a2
d_gt;
d = pdist([a2_computed;a1_computed]);
error = d - d_gt;

%% Initial error setup
figure(1)

image_error = deg2rad(2.367);%Set 3 Image Error
%image_error = deg2rad(5.586);%Set 2 Image Error
% image_error = deg2rad(5.586);%Set 6 Image Error

o = [0 0 0];
scatter3(o(1),o(2),o(3));
qiver_T = quiver3(o(1),o(2),o(3),T(1),T(2),T(3),0);
qiver_T.Color = 'green';
qiver_a1 = quiver3(T(1),T(2),T(3),a1(1),a1(2),a1(3),0);
qiver_a1.Color = 'blue';
qiver_a2 = quiver3(T(1),T(2),T(3),a2(1),a2(2),a2(3),0);
qiver_a2.Color = 'cyan';
qiver_b1 = quiver3(o(1),o(2),o(3),b1(1),b1(2),b1(3),0);
qiver_b1.Color = [0.5,0.5,0.5];
qiver_b2 = quiver3(o(1),o(2),o(3),b2(1),b2(2),b2(3),0);
qiver_b2.Color = 'red';

%% 1. Translation error volumes

%Define the number of points in the error volume/cloud
azStepSize = deg2rad(360)/numStepsAz;
elStepSize = deg2rad(180)/numStepsEl;
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

%% 2. Create distance error volume
if(compute_error_cloud)
tic
coneStepSize = deg2rad(360)/numSteps;
T_vec = [T(1); T(2); T(3)];
dL_error_array = [];
L_a1_array = [];
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
    
    d_alpha1 = alpha1 - alpha1_dash;
    d_beta1 = beta1 - beta1_dash;
    d_gamma1 = gamma1 - gamma1_dash;
    dT_dash = norm(T)-norm(T_dash);
    
    %approximation of the change in L (see eqn11 of paper)
    dL_a1_approx = norm(T)*((d_alpha1/sin(gamma1))-(d_gamma1*sin(alpha1)/sin(gamma1)^2))+dT_dash*(sin(alpha1)/sin(gamma1));

    %using the sine rule to compute the length of a1 and a2
    a1_len_dash = sin(alpha1_dash)/sin(gamma1_dash) * norm(T_dash);
    a2_len_dash = sin(alpha2_dash)/sin(gamma2_dash) * norm(T_dash);    
    b1_len_dash = sin(beta1_dash)/sin(gamma1_dash) * norm(T_dash);
    b2_len_dash = sin(beta2_dash)/sin(gamma2_dash) * norm(T_dash); 
    
    %Add image cone error b1
    %a1
    a1_norm = a1/norm(a1);
    dL_a1 =  a1_len - a1_len_dash;
    L_a1_array(i) = a1_len;
    dL_error_array(i) = dL_a1-dL_a1_approx;
    dG_a1 = (cross(a1_norm,cross(T_dash_norm,a1_norm)))/norm(cross(a1_norm,cross(T_dash_norm,a1_norm)))*sin(image_error)*(L_a1+dL_a1);
    dG_img_error_a1 = a1_norm*(L_a1+dL_a1) + dG_a1;

    %a2
    a2_norm = a2/norm(a2);
    dL_a2 =  a2_len - a2_len_dash;
    dG_a2 = (cross(a2_norm,cross(T_dash_norm,a2_norm)))/norm(cross(a2_norm,cross(T_dash_norm,a2_norm)))*sin(image_error)*(L_a2+dL_a2);
    dG_img_error_a2 = a2_norm*(L_a2+dL_a2) + dG_a2;
    %b1
    b1_norm = b1/norm(b1);
    dL_b1 =  b1_len - b1_len_dash;
    dG_b1 = (cross(b1_norm,cross(T_dash_norm,b1_norm)))/norm(cross(b1_norm,cross(T_dash_norm,b1_norm)))*sin(image_error)*(L_b1+dL_b1);
    dG_img_error_b1 = b1_norm*(L_b1+dL_b1) + dG_b1;
    %b2
    b2_norm = b2/norm(b2);
    dL_b2 =  b2_len - b2_len_dash;
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

% 3. Point cloud plotting
if(plot_cloud_points)
figure(1)
mesh(Sx_1,Sy_1,Sz_1);
hold on;
mesh(Sx_2,Sy_2,Sz_2);
%a1
scatter3(error_vol_x_a1,error_vol_y_a1,error_vol_z_a1, 'MarkerEdgeColor','b','MarkerFaceColor', 'b');
%a2
scatter3(error_vol_x_a2,error_vol_y_a2,error_vol_z_a2, 'MarkerEdgeColor','cyan','MarkerFaceColor', 'cyan');
%b1
scatter3(error_vol_x_b1,error_vol_y_b1,error_vol_z_b1, 'MarkerEdgeColor',[0.5,0.5,0.5],'MarkerFaceColor', [0.5,0.5,0.5]);
%b2
scatter3(error_vol_x_b2,error_vol_y_b2,error_vol_z_b2, 'MarkerEdgeColor','black','MarkerFaceColor', 'red');

axis equal
hold off;
end
% fprintf('Point cloud time: ')
PtCT=toc;
end
%%
%Compute the max and minimum distances from the point clouds.
if(compute_cloud_range)
tic
error_a1 = cat(2,error_vol_x_a1,error_vol_y_a1,error_vol_z_a1);
error_a2 = cat(2,error_vol_x_a2,error_vol_y_a2,error_vol_z_a2);
error_b1 = cat(2,error_vol_x_b1,error_vol_y_b1,error_vol_z_b1);
error_b2 = cat(2,error_vol_x_b2,error_vol_y_b2,error_vol_z_b2);

error_ab1 = [error_a1;error_b1];
error_ab2 = [error_a2;error_b2];
d_max_cloud = 0;
d_min_cloud = 0;
first = true;

reverseStr = '';

for i = 1:length(error_ab1)
    for j = 1:length(error_ab2)
        cloud_point_dist = pdist([error_ab1(i, :);error_ab2(j, :)]);
        if(first == true)
            d_max_cloud = cloud_point_dist;
            d_min_cloud = cloud_point_dist;
            first = false;
        else
            if cloud_point_dist > d_max_cloud
                d_max_cloud = cloud_point_dist;
            end
            if cloud_point_dist < d_min_cloud
                d_min_cloud = cloud_point_dist;
            end
        end
               
    end
    
    % Display the progress
   percentDone = 100 * i / length(error_ab1);
   msg = sprintf('Cloud range computation percent done: %3.1f', percentDone); %Don't forget this semicolon
   fprintf([reverseStr, msg]);
   reverseStr = repmat(sprintf('\b'), 1, length(msg));
end

range_point_cloud = [d_min_cloud, d_max_cloud];
NSd_C = d_gt\(d_max_cloud-d_min_cloud);
% fprintf('Point cloud min/max distance calculation time: ')
PtCMM_T = toc;


%Plot the percentage error of the approximation of L
figure(2)
plot(dL_error_array./L_a1_array*100.0)
xlabel('Error Cloud Pairs', 'FontSize', 14);
ylabel('L error (%)', 'FontSize', 14);

end

%% Compute the sphere approximation
if(compute_sphere_approx)
tic

dL_a1_max = 0;
dL_a2_max = 0;
dL_b1_max = 0;
dL_b2_max = 0;
for i=1:length(T_dash_X)
    %Reformate the tranlation vectors
    T_dash = [T_dash_X(i); T_dash_Y(i);T_dash_Z(i)];
    T_dash_norm = T_dash/norm(T_dash);
    %Compute dAlpha, dBeta, dGama, dL
    %Compute the new angles due to errors
    alpha1_dash = acos(dot(T_dash_norm,b1/norm(b1))); %Alpa_n calc from page 3 in paper
    alpha2_dash = acos(dot(T_dash_norm,b2/norm(b2)));
    beta1_dash = acos(-dot(T_dash_norm,a1/norm(a1)));
    beta2_dash = acos(-dot(T_dash_norm,a2/norm(a2)));
    
    %gamma will be useful for the sine rule in the next step
    gamma1_dash = pi()-alpha1_dash-beta1_dash;
    gamma2_dash = pi()-alpha2_dash-beta2_dash;

    %using the sine rule to compute the length of a1 and a2
    a1_len_dash = sin(alpha1_dash)/sin(gamma1_dash) * norm(T_dash); %eq 2 in paper
    dL_a1 =  a1_len - a1_len_dash;
    if(dL_a1>dL_a1_max)
        dL_a1_max = dL_a1;
    end

    a2_len_dash = sin(alpha2_dash)/sin(gamma2_dash) * norm(T_dash); 
    dL_a2 =  a2_len - a2_len_dash;
    if(dL_a2>dL_a2_max)
        dL_a2_max = dL_a2;
    end
    
    b1_len_dash = sin(beta1_dash)/sin(gamma1_dash) * norm(T_dash);
    dL_b1 =  b1_len - b1_len_dash;
    if(dL_b1>dL_b1_max)
        dL_b1_max = dL_b1;
    end
    
    b2_len_dash = sin(beta2_dash)/sin(gamma2_dash) * norm(T_dash); 
    dL_b2 =  b2_len - b2_len_dash;
    if(dL_b2>dL_b2_max)
        dL_b2_max = dL_b2;
    end
end
theta_T = atan(2*dT/norm(T)); %eq 6 in code

ra1_T = (L_a1+dL_a1_max)*sin(theta_T);
ra2_T = (L_a2+dL_a2_max)*sin(theta_T);
rb1_T = (L_b1+dL_b1_max)*sin(theta_T);
rb2_T = (L_b2+dL_b2_max)*sin(theta_T);

ra1_psi = (L_a1+dL_a1_max)*sin(image_error);
ra2_psi = (L_a2+dL_a2_max)*sin(image_error);
rb1_psi = (L_b1+dL_b1_max)*sin(image_error);
rb2_psi = (L_b2+dL_b2_max)*sin(image_error);

R_a1 = ra1_T + ra1_psi;
R_a2 = ra2_T + ra2_psi;
R_b1 = rb1_T + rb1_psi;
R_b2 = rb2_T + rb2_psi;

[x_a1,y_a1,z_a1] = sphere_points(R_a1,T(1) + a1_computed(1),T(2) + a1_computed(2),T(3) + a1_computed(3));
[x_a2,y_a2,z_a2] = sphere_points(R_a2,T(1) + a2_computed(1),T(2) + a2_computed(2),T(3) + a2_computed(3));
[x_b1,y_b1,z_b1] = sphere_points(R_b1,b1_computed(1),b1_computed(2),b1_computed(3));
[x_b2,y_b2,z_b2] = sphere_points(R_b2,b2_computed(1),b2_computed(2),b2_computed(3));

% fprintf('Approximation cloud calc time :')
PtAT=toc;
tic
d_a_min = norm(a2_computed-a1_computed)-R_a1-R_a2;
d_a_max = norm(a2_computed-a1_computed)+R_a1+R_a2;
d_b_min = norm(b2_computed-b1_computed)-R_b1-R_b2;
d_b_max = norm(b2_computed-b1_computed)+R_b1+R_b2;

d_min_sphere_approx = max([d_a_min,d_b_min]);
d_max_sphere_approx = min([d_a_max,d_b_max]);
range_sphere_approx = [d_min_sphere_approx, d_max_sphere_approx];
NSd_A = d_gt\(d_max_sphere_approx-d_min_sphere_approx);
% fprintf('Approximation cloud min/max calc time: ')
PtAMMT=toc;
% d_gt
% numStepsEl
% numStepsAz
% numSteps 
% T_Angle
% T_dist
% dT 
% Gap_size
%Plot spheres
figure(1);
view(-29,20)
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
%% Run Details table
    PointC_Info={'TAngle';'Gap';'Translation';'E1Steps';'AzSteps';'Steps';'ImageError';'dT'};
    Point_Cloud_Info = table(T_Angle,Gap_size,T_dist,numStepsEl,numStepsAz,numSteps,image_error,dT,'VariableNames',PointC_Info)
%% Point Cloud Comparison table
    PointC_Comparison={'PtCTime';'PtCMinMaxTime';'PtATime';'PtAMinMaxTime';'PtCMin';'PtCMax';'PtAMin';'PtAMax';'PtCRange';'PtARange';'RangeDiffs';'NSd_C';'NSd_A'};
    Point_Cloud_Data = table(PtCT,PtCMM_T,PtAT,PtAMMT,d_min_cloud,d_max_cloud,d_min_sphere_approx,d_max_sphere_approx,d_max_cloud-d_min_cloud,d_max_sphere_approx-d_min_sphere_approx,-(d_max_cloud-d_min_cloud)+(d_max_sphere_approx-d_min_sphere_approx),NSd_C,NSd_A,'VariableNames',PointC_Comparison)
end 




