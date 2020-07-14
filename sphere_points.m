function [x,y,z] = sphere_points(r,xc,yc,zc)
    N = 20;
    thetavec = linspace(0,pi,N);
    phivec = linspace(0,2*pi,2*N);
    [th, ph] = meshgrid(thetavec,phivec);
    R = r*ones(size(th)); % should be your R(theta,phi) surface in general

    x = xc + R.*sin(th).*cos(ph);
    y = yc + R.*sin(th).*sin(ph);
    z = zc + R.*cos(th);
end