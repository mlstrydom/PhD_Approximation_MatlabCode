clear all;
close all;
plot_figs = 0;
T_dist = 8;%Translation Distance
%Noise values
%---------------------------Inputs ---------------------------------
dT = 0.0698; %Translation Error
image_error = deg2rad(5.586);
%------------------------------------------------------------------
error_a_array(1) = 0;
error_b_array(1) = 0;
error_array(1) = 0;

%%
for dir = deg2rad(0):deg2rad(15):deg2rad(90)
for run_num = 1:2

    %Athroscope positions (current = A2, previous = A1).
    %Assumes that the rotation has been accounted for (i.e. pure translation).
    
    A1_tmp = [-T_dist*sin(dir) -T_dist*cos(dir) 0];
    A1 = [-T_dist*sin(dir) -T_dist*cos(dir) 0]-A1_tmp;
    A2 = [0 0 0]-A1_tmp;

    % Select two points to compute the distance (e.g. two points detected on the knee gap).
    G1 = [-1 3 0];
    G2 = [1 2 1];
    G_dir = (G2-G1)/norm(G2-G1);
    G_dist = 4;
    G2 = G1+G_dist*G_dir;
    G1 = G1-A1_tmp;
    G2 = G2-A1_tmp;

    %plot GT positions
    if plot_figs==1
        figure(1)
        scatter3(A1(1), A1(2), A1(3), 'b','LineWidth',3)
        hold on;
        scatter3(A2(1), A2(2), A2(3), 'b','LineWidth',3)
        scatter3(G1(1), G1(2), G1(3), 'r','LineWidth',3)
        scatter3(G2(1), G2(2), G2(3), 'r','LineWidth',3)
        A = [A1;A2];
        G = [G1;G2];
        plot3(A(:,1),A(:,2),A(:,3),'b','LineWidth',3)
        plot3(G(:,1),G(:,2),G(:,3),'r','LineWidth',3)
        axis equal
        xyzlabel()
        set(gca,'fontsize',20)
    end

    %Used to add translation noise
    noise_T = [-dT dT];
    n=1;

    %Compute the ground truth distance
    d_gt = pdist([G1;G2]);
    %Translation of athroscope tip
    T = A2 - A1;
    T = T + rand(n,1)*range(noise_T)+min(noise_T);

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

    %Used to add translation noise

    noise_image = [deg2rad(-image_error) deg2rad(image_error)];

    %Compute the angles we know
    alpha1 = acos(dot(T/norm(T),b1/norm(b1)));
    alpha1 = alpha1 + rand(n,1)*range(noise_image)+min(noise_image);
    alpha2 = acos(dot(T/norm(T),b2/norm(b2)));
    alpha2 = alpha2 + rand(n,1)*range(noise_image)+min(noise_image);
    beta1 = acos(dot(-T/norm(T),a1/norm(a1)));
    beta1 = beta1 + rand(n,1)*range(noise_image)+min(noise_image);
    beta2 = acos(dot(-T/norm(T),a2/norm(a2)));
    beta2 = beta2 + rand(n,1)*range(noise_image)+min(noise_image);

    %gamma will be useful for the sine rule in the next step
    gamma1 = pi()-alpha1-beta1;
    gamma2 = pi()-alpha2-beta2;

    %using the sine rule to compute the length of a1 and a2
    a1_len = sin(alpha1)/sin(gamma1) * norm(T);
    a2_len = sin(alpha2)/sin(gamma2) * norm(T);
    b1_len = sin(beta1)/sin(gamma1) * norm(T);
    b2_len = sin(beta2)/sin(gamma2) * norm(T);

    %now compute the a1 and a2 vectors
    a1_computed = a1_len * (a1/norm(a1));
    a2_computed = a2_len * (a2/norm(a2));
    b1_computed = b1_len * (b1/norm(b1));
    b2_computed = b2_len * (b2/norm(b2));

    %Compute the distance between a1 and a2
    d_gt
    d_a = pdist([a2_computed;a1_computed]);
    d_b = pdist([b1_computed;b2_computed]);
    error_a = d_a - d_gt;
    error_b = d_b - d_gt;

    % Initial error setup
    if plot_figs==1
        figure(1)
        o = [0 0 0];
        scatter3(o(1),o(2),o(3))
        quiver3(A1(1),A1(2),A1(3),A2(1),A2(2),A2(3),0,'LineWidth',3)
        quiver3(A2(1),A2(2),A2(3),a1(1),a1(2),a1(3),0,'LineWidth',3)
        quiver3(A2(1),A2(2),A2(3),a2(1),a2(2),a2(3),0,'LineWidth',3)
        quiver3(A1(1),A1(2),A1(3),b1(1),b1(2),b1(3),0,'LineWidth',3)
        quiver3(A1(1),A1(2),A1(3),b2(1),b2(2),b2(3),0,'LineWidth',3)
    end

    % Translation error volumes
    %Define the number of points in the error volume/cloud
    numStepsEl = 8;
    numStepsAz = 16;
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
    T_dash_X = SX_2(i1)-SX_1(i2);
    T_dash_Y = SY_2(i1)-SY_1(i2);
    T_dash_Z = SZ_2(i1)-SZ_1(i2);
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

    % Create distance error volume
    numSteps = 4;
    coneStepSize = deg2rad(360)/numSteps;
    T_vec = [T(1); T(2); T(3)];

    for i=1:length(T_dash_X)
        %Reformate the tranlation vectors
        T_dash = [T_dash_X(i); T_dash_Y(i);T_dash_Z(i)];
        S1 = [SX_1_array(i); SY_1_array(i);SZ_1_array(i)];
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
        dG_img_error_a1 = a1_norm*(L_a1+dL_a1) + dG_a1 + T;
        %a2
        a2_norm = a2/norm(a2);
        dL_a2 =  a2_len - a2_len_dash;
        L_a2 = a2_len;
        dG_a2 = (cross(a2_norm,cross(T_dash_norm,a2_norm)))/norm(cross(a2_norm,cross(T_dash_norm,a2_norm)))*sin(image_error)*(L_a2+dL_a2);
        dG_img_error_a2 = a2_norm*(L_a2+dL_a2) + dG_a2 + T;
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
            v_rot_a1 = transpose(S1)+rodrigues_rot(dG_img_error_a1,a1_norm,theta);
            v_rot_x_a1(i,j) = v_rot_a1(1);
            v_rot_y_a1(i,j) = v_rot_a1(2);
            v_rot_z_a1(i,j) = v_rot_a1(3);
            %a2
            v_rot_a2 = transpose(S1)+rodrigues_rot(dG_img_error_a2,a2_norm,theta);
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

    %
    %Compute noise to signal of gap distance
    pairWiseDistances_a = pdist2([error_vol_x_a1, error_vol_y_a1, error_vol_z_a1], [error_vol_x_a2, error_vol_y_a2,error_vol_z_a2]);
    [maxDistance_a, maxLinearIndex_a] = max(pairWiseDistances_a(:));
    [minDistance_a, minLinearIndex_a] = min(pairWiseDistances_a(:));

    %noiseToSignal_a(run_num) = (maxDistance_a-minDistance_a)/d_a;
    noiseToSignal_a(run_num) = d_a/(maxDistance_a-minDistance_a);
    %SNR_a(run_num)=1/noiseToSignal_a(run_num);

    pairWiseDistances_b = pdist2([error_vol_x_b1, error_vol_y_b1, error_vol_z_b1], [error_vol_x_b2, error_vol_y_b2,error_vol_z_b2]);
    [maxDistance_b, maxLinearIndex_b] = max(pairWiseDistances_b(:));
    [minDistance_b, minLinearIndex_b] = min(pairWiseDistances_b(:));
    %noiseToSignal_b(run_num) = (maxDistance_b-minDistance_b)/d_b;
    noiseToSignal_b(run_num) = d_b/(maxDistance_b-minDistance_b);
    %SNR_b(run_num)=1/noiseToSignal_b(run_num);
    % Plotting
    if plot_figs==1
        figure(1)
        mesh(Sx_1,Sy_1,Sz_1)
        hold on;
        mesh(Sx_2,Sy_2,Sz_2)
        %a1
        scatter3(error_vol_x_a1,error_vol_y_a1,error_vol_z_a1)
        %a2
        scatter3(error_vol_x_a2,error_vol_y_a2,error_vol_z_a2)
        %b1
        scatter3(error_vol_x_b1,error_vol_y_b1,error_vol_z_b1)
        %b2
        scatter3(error_vol_x_b2,error_vol_y_b2,error_vol_z_b2)
        % DT = delaunayTriangulation(x,y,z);  %# Create the tetrahedral mesh
        % tetramesh(DT);

        hold off;
        
    end
    %Arrays
    error_a_array(run_num) = error_a;
    error_b_array(run_num) = error_b;
    error_array(run_num) = min(error_a,error_b);    
end
display('Angle run: ')
rad2deg(dir)

mean_error = mean(error_array)
std_error = std(error_array)

end
%%
figure()
histogram(error_array)


%% 
%Translation impact
clear all;
plot_figs = 0;
run_num = 1;
for T_dist = 0.1:0.1:8 %Vary the translation distance from 0.5 to 8mm
    T_dist
    %Noise values
     dT = 0.0698;
    %dT = 0.0698;
    image_error = deg2rad(5.586);
%     image_error = deg2rad(5.586);
    dir = deg2rad(45);
    %Athroscope positions (current = A2, previous = A1).
    %Assumes that the rotation has been accounted for (i.e. pure translation).
    
    A1_tmp = [-T_dist*sin(dir) -T_dist*cos(dir) 0];
    A1 = [-T_dist*sin(dir) -T_dist*cos(dir) 0]-A1_tmp;
    A2 = [0 0 0]-A1_tmp;

    %Select two points to compute the distance (e.g. two points detected on the knee gap).
    G1 = [-1 3 0];
    G2 = [1 2 1];
    G_dir = (G2-G1)/norm(G2-G1);
    G_dist = 4;
    G2 = G1+G_dist*G_dir;
    
    G1 = G1-A1_tmp;
    G2 = G2-A1_tmp;

    %plot GT positions
    if plot_figs==1
        figure(1)
        scatter3(A1(1), A1(2), A1(3), 'b','LineWidth',3)
        hold on;
        scatter3(A2(1), A2(2), A2(3), 'b','LineWidth',3)
        scatter3(G1(1), G1(2), G1(3), 'r','LineWidth',3)
        scatter3(G2(1), G2(2), G2(3), 'r','LineWidth',3)
        A = [A1;A2];
        G = [G1;G2];
        plot3(A(:,1),A(:,2),A(:,3),'b','LineWidth',3)
        plot3(G(:,1),G(:,2),G(:,3),'r','LineWidth',3)
        axis([0 8 0 15])
        %axis equal
        xyzlabel()
        set(gca,'fontsize',20)
    end

    %Used to add translation noise
    noise_T = [-dT dT];
    n=1;

    %Compute the ground truth distance
    d_gt = pdist([G1;G2]);
    %Translation of athroscope tip
    T = A2 - A1;
    T = T + rand(n,1)*range(noise_T)+min(noise_T);

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

    %Used to add translation noise

    noise_image = [deg2rad(-image_error) deg2rad(image_error)];

    %Compute the angles we know
    alpha1 = acos(dot(T/norm(T),b1/norm(b1)));
    alpha1 = alpha1 + rand(n,1)*range(noise_image)+min(noise_image);%Add rotational and image noise 
    alpha2 = acos(dot(T/norm(T),b2/norm(b2)));
    alpha2 = alpha2 + rand(n,1)*range(noise_image)+min(noise_image);
    beta1 = acos(dot(-T/norm(T),a1/norm(a1)));
    beta1 = beta1 + rand(n,1)*range(noise_image)+min(noise_image);
    beta2 = acos(dot(-T/norm(T),a2/norm(a2)));
    beta2 = beta2 + rand(n,1)*range(noise_image)+min(noise_image);

    %gamma will be useful for the sine rule in the next step
    gamma1 = pi()-alpha1-beta1;
    gamma2 = pi()-alpha2-beta2;

    %using the sine rule to compute the length of a1 and a2
    a1_len = sin(alpha1)/sin(gamma1) * norm(T);
    a2_len = sin(alpha2)/sin(gamma2) * norm(T);
    b1_len = sin(beta1)/sin(gamma1) * norm(T);
    b2_len = sin(beta2)/sin(gamma2) * norm(T);

    %now compute the a1 and a2 vectors
    a1_computed = a1_len * (a1/norm(a1));
    a2_computed = a2_len * (a2/norm(a2));
    b1_computed = b1_len * (b1/norm(b1));
    b2_computed = b2_len * (b2/norm(b2));

    %Compute the distance between a1 and a2
    d_gt;
    d_a = pdist([a2_computed;a1_computed]);
    d_b = pdist([b1_computed;b2_computed]);
    error_a = d_a - d_gt;
    error_b = d_b - d_gt;

    % Initial error setup
    if plot_figs==1
        figure(1)
        o = [0 0 0];
        scatter3(o(1),o(2),o(3))
        quiver3(A1(1),A1(2),A1(3),A2(1),A2(2),A2(3),0,'LineWidth',3)
        quiver3(A2(1),A2(2),A2(3),a1(1),a1(2),a1(3),0,'LineWidth',3)
        quiver3(A2(1),A2(2),A2(3),a2(1),a2(2),a2(3),0,'LineWidth',3)
        quiver3(A1(1),A1(2),A1(3),b1(1),b1(2),b1(3),0,'LineWidth',3)
        quiver3(A1(1),A1(2),A1(3),b2(1),b2(2),b2(3),0,'LineWidth',3)
    end

    % Translation error volumes
    %Define the number of points in the error volume/cloud
    numStepsEl = 8;
    numStepsAz = 16;
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
    T_dash_X = SX_2(i1)-SX_1(i2);
    T_dash_Y = SY_2(i1)-SY_1(i2);
    T_dash_Z = SZ_2(i1)-SZ_1(i2);
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

    % Create distance error volume
    numSteps = 4;
    coneStepSize = deg2rad(360)/numSteps;
    T_vec = [T(1); T(2); T(3)];

    for i=1:length(T_dash_X)
        %Reformate the tranlation vectors
        T_dash = [T_dash_X(i); T_dash_Y(i);T_dash_Z(i)];
        S1 = [SX_1_array(i); SY_1_array(i);SZ_1_array(i)];
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
        dG_img_error_a1 = a1_norm*(L_a1+dL_a1) + dG_a1 + T;
        %a2
        a2_norm = a2/norm(a2);
        dL_a2 =  a2_len - a2_len_dash;
        L_a2 = a2_len;
        dG_a2 = (cross(a2_norm,cross(T_dash_norm,a2_norm)))/norm(cross(a2_norm,cross(T_dash_norm,a2_norm)))*sin(image_error)*(L_a2+dL_a2);
        dG_img_error_a2 = a2_norm*(L_a2+dL_a2) + dG_a2 + T;
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
            v_rot_a1 = transpose(S1)+rodrigues_rot(dG_img_error_a1,a1_norm,theta);
            v_rot_x_a1(i,j) = v_rot_a1(1);
            v_rot_y_a1(i,j) = v_rot_a1(2);
            v_rot_z_a1(i,j) = v_rot_a1(3);
            %a2
            v_rot_a2 = transpose(S1)+rodrigues_rot(dG_img_error_a2,a2_norm,theta);
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

    %
    %Compute noise to signal of gap distance
    pairWiseDistances_a = pdist2([error_vol_x_a1, error_vol_y_a1, error_vol_z_a1], [error_vol_x_a2, error_vol_y_a2,error_vol_z_a2]);
    [maxDistance_a, maxLinearIndex_a] = max(pairWiseDistances_a(:));
    [minDistance_a, minLinearIndex_a] = min(pairWiseDistances_a(:));

%     noiseToSignal_a(run_num) = (maxDistance_a-minDistance_a)/d_a;
    noiseToSignal_a(run_num) = d_a/(maxDistance_a-minDistance_a);

    pairWiseDistances_b = pdist2([error_vol_x_b1, error_vol_y_b1, error_vol_z_b1], [error_vol_x_b2, error_vol_y_b2,error_vol_z_b2]);
    [maxDistance_b, maxLinearIndex_b] = max(pairWiseDistances_b(:));
    [minDistance_b, minLinearIndex_b] = min(pairWiseDistances_b(:));
%     noiseToSignal_b(run_num) = (maxDistance_b-minDistance_b)/d_b;
    noiseToSignal_b(run_num) = d_b/(maxDistance_b-minDistance_b);

    % Plotting
    if plot_figs==1
        figure(1)
        mesh(Sx_1,Sy_1,Sz_1)
        hold on;
        mesh(Sx_2,Sy_2,Sz_2)
        %a1
        scatter3(error_vol_x_a1,error_vol_y_a1,error_vol_z_a1)
        %a2
        scatter3(error_vol_x_a2,error_vol_y_a2,error_vol_z_a2)
        %b1
        scatter3(error_vol_x_b1,error_vol_y_b1,error_vol_z_b1)
        %b2
        scatter3(error_vol_x_b2,error_vol_y_b2,error_vol_z_b2)
        % DT = delaunayTriangulation(x,y,z);  %# Create the tetrahedral mesh
        % tetramesh(DT);

        hold off;

        %Arrays
        error_array(run_num) = min(error_a,error_b);
    end
    run_num = run_num + 1;
end

%%
run_num

T_dist = 0.1:0.1:8;
size(T_dist)
size(noiseToSignal_a)
figure(3)
plot(T_dist, noiseToSignal_a,'b','LineWidth',3)
hold on;
plot(T_dist, noiseToSignal_b,'g','LineWidth',3)
xlabel('Translation Distance (mm)')
ylabel('Signal to Noise Ratio')
        axis([0 8 0 15])
legend('SNR(b_n)','SNR(a_n)')
set(gca,'fontsize',14)


