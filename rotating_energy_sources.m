function MVTS=rotating_energy_sources(T, options)

%% Initializing the parameters. 
opts = struct('N', 20, 'centers', [ [15; 15], [5; 5] ], ...
	'angularVelocities', [2*pi/20, 2*pi/5], ...
	'radii', [5, 5], 'sourceEnergy', 10, 'bandwidth', 25, ...
	'phases', [2*pi, 2*pi].*rand(1,2), ...
	'noiseVar', 1, 'noiseMean', 0, 'noiseCorrCoeff', 0);

if nargin > 1
	fields = fieldnames(options);
	for f=1:length(fields),
		opts.(fields{f}) = options.(fields{f});
	end
end

disp(opts);

N = opts.N;

%% The centers are on the columns.
centers = opts.centers;
angularVelocities = opts.angularVelocities;
radii = opts.radii;

sourceEnergy = opts.sourceEnergy;
bandwidth    = opts.bandwidth;

phases = opts.phases; 


sourcePositions = zeros(size(centers)); 
%% This is the variable that the function returns.
MVTS=zeros(N, N, T); 


noiseVar = opts.noiseVar;
noiseMean = opts.noiseMean;
noiseCorrCoeff = opts.noiseCorrCoeff;

if( noiseCorrCoeff > 0 )
	Wprev=randn(N,N);
end

for t=1:T,
	phases = mod(phases + angularVelocities, 2*pi);
	for idx=1:size(sourcePositions, 2),
		sourcePositions(:, idx) = centers(:, idx) + radii(idx)*[cos(phases(idx)); sin(phases(idx))];
	end

	spatial = zeros(N, N);
	
	for idx=1:size(sourcePositions, 2),

		mu    = sourcePositions(:, idx);
		Sigma = (bandwidth/2)*eye(2);

		[X,Y] = meshgrid(0:N-1, 0:N-1);
		C  = [X(:) Y(:)];
		downLeft = reshape( mvncdf(C  , mu', Sigma), N, N);
		upRight  = reshape( mvncdf(C+1, mu', Sigma), N, N);
		C  = [(X(:)+1) Y(:)];
		downRight = reshape( mvncdf(C, mu', Sigma), N, N);
		C  = [X(:) (Y(:)+1)];
		upLeft    = reshape( mvncdf(C, mu', Sigma), N, N);

		dynamic = sourceEnergy*(upRight+downLeft-downRight-upLeft)*(2*pi*bandwidth);
		spatial = spatial + dynamic;
	end

	%%Add noise.
	W = randn(N, N);
	if( noiseCorrCoeff > 0 )
		red     = noiseCorrCoeff.*Wprev + sqrt(1-noiseCorrCoeff^2).*W;
		spatial = spatial + sqrt(noiseVar).*red + noiseMean;
	else	
		spatial = spatial + sqrt(noiseVar).*W + noiseMean;
	end

	MVTS(:,:,t) = spatial;
end

MVTS = MVTS - repmat(mean(MVTS, 3), [1, 1, T]);

%for t=1:T,
%	minVal  = prctile(MVTS(:), 1);
%	maxVal  = prctile(MVTS(:), 99);
%	spatial = MVTS(:, :, t);
%	imagesc(spatial, [minVal, maxVal]);
%	colorbar;
%	pause(.3);
%
%end

end
