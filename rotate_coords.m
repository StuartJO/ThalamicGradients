function newcoords = rotate_coords(coords,azel,alpha,origin)
%ROTATE Rotate objects about specified origin and direction.
%   ROTATE(H,[THETA PHI],ALPHA) rotates the objects with handles H
%   through angle ALPHA about an axis described by the 2-element
%   direction vector [THETA PHI] (spherical coordinates).  
%   All the angles are in degrees.  The handles in H must be children
%   of the same axes.
%
%   THETA is the angle in the xy plane counterclockwise from the
%   positive x axis.  PHI is the elevation of the direction vector
%   from the xy plane (see also SPH2CART).  Positive ALPHA is defined
%   as the righthand-rule angle about the direction vector as it
%   extends from the origin.
%
%   ROTATE(H,[X Y Z],ALPHA) rotates the objects about the direction
%   vector [X Y Z] (Cartesian coordinates). The direction vector
%   is the vector from the center of the plot box to (X,Y,Z).
%
%   ROTATE(...,ORIGIN) uses the point ORIGIN = [x0,y0,z0] as
%   the center of rotation instead of the center of the plot box.
%
%   See also SPH2CART, CART2SPH.

%   Copyright 1984-2017 The MathWorks, Inc.

% Determine the default origin (center of plot box).

% find unit vector for axis of rotation
if numel(azel) == 2 % theta, phi
    theta = pi*azel(1)/180;
    phi = pi*azel(2)/180;
    u = [cos(phi)*cos(theta); cos(phi)*sin(theta); sin(phi)];
elseif numel(azel) == 3 % direction vector
    u = azel(:)/norm(azel);
end

alph = alpha*pi/180;
cosa = cos(alph);
sina = sin(alph);
vera = 1 - cosa;
x = u(1);
y = u(2);
z = u(3);
rot = [cosa+x^2*vera x*y*vera-z*sina x*z*vera+y*sina; ...
       x*y*vera+z*sina cosa+y^2*vera y*z*vera-x*sina; ...
       x*z*vera-y*sina y*z*vera+x*sina cosa+z^2*vera]';


  skip = 0;

       x = coords(:,1); y = coords(:,2); 
       if size(coords,2)>2
          z = coords(:,3);
       else
          z = [];
       end
       

    if isempty(z)
       z = -origin(3)*ones(size(y));
    end
    [m,n] = size(z);
    if numel(x) < m*n
      [x,y] = meshgrid(x,y);
    end
  
  if ~skip
    [m,n] = size(x);
    newxyz = [x(:)-origin(1), y(:)-origin(2), z(:)-origin(3)];
    newxyz = newxyz*rot;
    newx = origin(1) + reshape(newxyz(:,1),m,n);
    newy = origin(2) + reshape(newxyz(:,2),m,n);
    newz = origin(3) + reshape(newxyz(:,3),m,n);

    newcoords = [newx,newy,newz];
  end



