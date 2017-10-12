if ~exist('model','var')
    model = GOLDModel('/Volumes/RadiativeTr/gold/');
end

lat0 = 40;
lon0 = -150;
maxT = 10*86400;
increment = 3600;
time = (0:3600:maxT)';

% Seems that an absolute tolerance of 1e-4 is appropriate, this corresponds
% to 0.0001 degrees, or 11 meters This also seems true of relative
% tolerance, at least for my test particle at -150,40

tic

func = @(time, position) model.PositionFlux(time, position);
options = odeset('RelTol',1e-4,'AbsTol',1e-4);
[T, X] = ode23(func,time,[lon0,lat0], options);
figure, plot( X(:,1),X(:,2) )
xlabel('longitude (degrees)')
ylabel('latitude (degrees)')

toc