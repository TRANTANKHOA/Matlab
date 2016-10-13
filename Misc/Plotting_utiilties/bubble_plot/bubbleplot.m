function h=bubbleplot(x,y,z,color,sf)
%function h=bubbleplot(x,y,z,color,sf)
%PURPOSE
%  Scatter plot in x and y using arbitrary-color "." as
%  plot symbol - whose point-by-point size is proportional
%  the the magnitude of z
%INPUT
%  x     - n-dimension vector
%  y     - n-dimension vector
%  z     - n-dimension vector (used to size the plot symbols)
%  color - plot symbol color (must be a 3-element vector
%          with elements in range 0 1) ................... default: alternating 10 colors
%  sf    - plot symbol size scale factor ................. default: 35 

n = nargin;
if n<5 | isempty(sf), sf = 20; end
if n<4 | isempty(color),
   myco=[ 0     0     1.00
          0     0.50  0
          1.00  0     0
          0     0.75  0.75
          0.75  0.75  0.75
          0.75  0     0.75
          0.25  1.00  0.25
          0.75  0.75  0
          0.25  0.25  0.25
          0.50  0.50  0.50 ];
else   
   col = color(:)';
   myco = [col;col;col;col;col;col;col;col;col;col];
end
scf=sf./max(z(:));
zp=round( z .*scf ); 
I = find(isnan(zp));
zp(I) = eps; 
for i = 1:length(x),
   cc = i;
   if cc>10,
      cc = mod(cc,10)+1;
   end
   plot(x(i),y(i),'.','MarkerSize',zp(i),'Color','blue');% myco(cc,:)
   hold on;
end
hold off
axis('tight')
ax=axis ;             % L R B T
xr=.03*(ax(2)-ax(1));
yr=.03*(ax(4)-ax(3));% scale axis based on data range
axis([ax(1)-xr ax(2)+2*xr ax(3)-yr ax(4)+2*yr]);
hold off;
h = gca;