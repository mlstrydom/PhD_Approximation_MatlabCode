%%

%Assumes that the rotation has been accounted for (i.e. pure translation).
% Zero degrees
% A1 = [0 0 0] %Athroscope positions (current = A2, previous = A1).
% A2 = [2 2 2]
% G1 = [0 3 0] %Select two points to compute the distance (e.g. two points detected on the knee gap).
% G2 = [2 5 2]
% %45 Degrees
% A1 = [0 0 0] %Athroscope positions (current = A2, previous = A1).
% A2 = [2 0 2]
% G1 = [0 0 0] %Select two points to compute the distance (e.g. two points detected on the knee gap).
% G2 = [2 0 0]

% % Old 45 Degrees (137Deg)
% A1 = [0 0 0] %Athroscope positions (current = A2, previous = A1).
% A2 = [-2 2 -0.7]
% G1 = [-1 3 0] %Select two points to compute the distance (e.g. two points detected on the knee gap).
% G2 = [1 2 -1]
% 135 Degrees
% A1 = [0 0 0] %Athroscope positions (current = A2, previous = A1).
% A2 = [-2 2 -0.87297]
% G1 = [-1 3 0] %Select two points to compute the distance (e.g. two points detected on the knee gap).
% G2 = [1 2 -1]
%90 Degrees
% A1 = [0 0 0] %Athroscope positions (current = A2, previous = A1).
% A2 = [2 2 0]
% G1 = [0 3 0] %Select two points to compute the distance (e.g. two points detected on the knee gap).
% G2 = [0 3 4]

At=A2-A1
Gt=G2-G1
AGDOT=dot(At,Gt)
AXG=cross(At,Gt)
AGXNorm=norm(AXG)
TransAngleCheck=atan2d(AGXNorm,AGDOT)
CosTAngel=acosd(AGDOT/(norm(At)*norm(Gt)))