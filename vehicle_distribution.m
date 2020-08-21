function [res] = vehicle_distribution (v)


%Road grid
% 
y1=0;y2=7;y3=426;y4=433;
 

x1=0;x2=7;x3=243;x4=250;
%coordinates in polar plan
xDelta=(x2-x1)+(x4-x3);
yDelta=(y2-y1)+(y4-y3);

%total area of the road grid
areaTotal=(433*250)-(419*236);

%spatial poisson process density
lambda = (1)/(2.5*v);
%spatial poisson process distribution 
numPoints = poissrnd(areaTotal*lambda);
xx=  xDelta*(rand(numPoints,1));  %angular coordinates
yy=  yDelta*(rand(numPoints,1));  %radial coordinates

%convert from polar to cartesian coordinates
[x,y]= pol2cart(xx,yy);

%Plotting
res = [x,y];
scatter(x,y);
xlabel('x');
ylabel('y');
axis square;



%figure;hold on;grid on;
%calculations for your first data
%plot(x1,y1,'-');
%calculations for your second data
%plot(x2,y2,'o');
%
end
