filename = 'ocean_geometry.nc';

lath = ncread(filename, 'lath');
lonh = ncread(filename, 'lonh');
wet = ncread(filename, 'wet');

lon =(-279:79)';
lat = (-77:85)';
[x0,y0] = ndgrid(lon,lat);
x0 = reshape(x0,numel(x0),1);
y0 = reshape(y0,numel(y0),1);

[X,Y] = ndgrid(lonh,lath);
wet_interp = griddedInterpolant(X,Y,wet);
iswet = wet_interp(x0,y0);
badindices = find(iswet == 0);
x0(badindices) = [];
y0(badindices) = [];


figure
pcolor(lonh,lath,wet')
shading flat

[x,y] = meshgrid(lonh,lath);
figure
scatter(reshape(x,1952*2880,1),reshape(y,1952*2880,1))