% Load the model, if necessary
if ~exist('model','var')
    model = GOLDModel('/Volumes/RadiativeTr/gold/');
end

% Now let's create a grid of ocean trajectories
filename = 'ocean_geometry.nc';

lath = ncread(filename, 'lath');
lonh = ncread(filename, 'lonh');
wet = ncread(filename, 'wet');

lon =(-279:10:79)';
lat = (-77:10:85)';
[x0,y0] = ndgrid(lon,lat);
x0 = reshape(x0,numel(x0),1);
y0 = reshape(y0,numel(y0),1);

[X,Y] = ndgrid(lonh,lath);
wet_interp = griddedInterpolant(X,Y,wet);
iswet = wet_interp(x0,y0);
badindices = find(iswet == 0);
x0(badindices) = [];
y0(badindices) = [];

initialConditions = [x0; y0];

% time
maxT = 100*86400;
increment = 3600;
time = (0:3600:maxT)';

tic
func = @(time, position) model.PositionFlux(time, position);
options = odeset('RelTol',1e-4,'AbsTol',1e-4);
[T, X] = ode23(func,time,initialConditions, options);
toc

lonIndices = (1:length(initialConditions)/2)';
latIndices = ((length(initialConditions)/2+1):length(initialConditions))';

lonT = X(:,lonIndices);
latT = X(:,latIndices);

figure, plot( lonT,latT )
xlabel('longitude (degrees)')
ylabel('latitude (degrees)')