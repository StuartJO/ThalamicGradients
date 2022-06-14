function [x,y] = rectangleNpoints(x_corner,y_corner,N,connect)

c(1,:) = [x_corner(1) y_corner(1)];
c(2,:) = [x_corner(1)+x_corner(2) y_corner(1)];
c(3,:) = [x_corner(1)+x_corner(2) y_corner(1)+y_corner(2)];
c(4,:) = [x_corner(1) y_corner(1)+y_corner(2)];

E{1} = [c(1,:);c(2,:)];
E{2} = [c(2,:);c(3,:)];
E{3} = [c(3,:);c(4,:)];
E{4} = [c(4,:);c(1,:)];

for i = 1:4
    m(i) = pdist(E{i});
end

Rm = m./sum(m);

NperSide = round(N.*Rm);

Rpoints = N-sum(NperSide);

x = [];
y = [];

for i = 1:4
    if i == 1
        extraPoints = Rpoints;
    else
        extraPoints = 0;
    end
    
    if i == 4
        ii = 1;
    else
        ii = i + 1;
    end
    xcoords = linspace(c(i,1),c(ii,1),NperSide(i)+1+extraPoints)';
    ycoords = linspace(c(i,2),c(ii,2),NperSide(i)+1+extraPoints)';
    
    if connect && i == 1
        xcoords = linspace(c(i,1),c(ii,1),NperSide(i)+extraPoints)';
        ycoords = linspace(c(i,2),c(ii,2),NperSide(i)+extraPoints)';
        
        x = [x; xcoords];
        y = [y; ycoords];
    else
        x = [x; xcoords(2:end)];
        y = [y; ycoords(2:end)];
    end
end