function [x,y] = circle(xCenter,yCenter,radius,nPoints)
angles = linspace(0, 2*pi, nPoints)'; % 720 is the total number of points
x = radius * cos(angles) + xCenter; 
y = radius * sin(angles) + yCenter;