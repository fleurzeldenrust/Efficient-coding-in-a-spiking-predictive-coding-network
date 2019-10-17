function p = pcolor_fleur(x,y,z)
% Same as pcolor, but than x and y are the middle of the columns/rows, not
% the edges


lx = length(x);
ly = length(y);

newx = zeros(1,lx+1);
newx(1) = x(1)-(x(2)-x(1))/2;
for n = 2:lx
    newx(n) = (x(n)+x(n-1))/2;
end
newx(lx+1) = x(lx)+(x(lx)-x(lx-1))/2;

newy = zeros(1,ly+1);
newy(1) = y(1)-(y(2)-y(1))/2;
for n = 2:ly
    newy(n) = (y(n)+y(n-1))/2;
end
newy(ly+1) = y(ly)+(y(ly)-y(ly-1))/2;

z = padarray(z,[1 1],'replicate','post');

p = pcolor(newx, newy, z);

