function [Z, trueZ] = rotating_energy_sources(T, options)

	%% Initializing the parameters to the default values. 
	opts = struct('N', 20, 'centers', [ [15; 15], [5; 5] ], ...
		'angularVelocities', [2*pi/20, 2*pi/5], ...
		'radii', [5, 5], 'sourceEnergy', 6.3, 'bandwidth', 50, ...
		'phases', [2*pi, 2*pi].*rand(1, 2), ...
		'noiseVar', 0, 'noiseMean', 0, 'noiseCorrCoeff', 0);

	if nargin > 1, %% Updating the parameters to the input values.
		fields = fieldnames(options);
		for f=1:length(fields),
			opts.(fields{f}) = options.(fields{f});
		end
	end

	%disp(opts);

	N = opts.N;

	%% The centers are on the columns.
	centers           = opts.centers;
	angularVelocities = opts.angularVelocities;
	radii             = opts.radii;
	phases            = opts.phases; 


	sourceEnergy = opts.sourceEnergy;
	bandwidth    = opts.bandwidth;


	sourcePositions = zeros(size(centers)); 
	nSources        = size(sourcePositions, 2);


	noiseVar       = opts.noiseVar;
	noiseMean      = opts.noiseMean;
	noiseCorrCoeff = opts.noiseCorrCoeff;

	%% Allocating the memory for the variables that the function returns.
	Z     = zeros(N, N, T); 
	trueZ = zeros(N, N, T, nSources); 


	if noiseCorrCoeff > 0,
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
			Sigma = (bandwidth/2);

			[X,Y] = meshgrid(0:N-1, 0:N-1);
			C  = [X(:) Y(:)];
			downLeft = reshape( mvncdf(C  , mu', Sigma*eye(2)), N, N);
			upRight  = reshape( mvncdf(C+1, mu', Sigma*eye(2)), N, N);
			C  = [(X(:)+1) Y(:)];
			downRight = reshape( mvncdf(C, mu', Sigma*eye(2)), N, N);
			C  = [X(:) (Y(:)+1)];
			upLeft    = reshape( mvncdf(C, mu', Sigma*eye(2)), N, N);

			dynamic = sourceEnergy*(upRight+downLeft-downRight-upLeft)*(2*pi*Sigma);
			spatial = spatial + dynamic;

			trueZ(:, :, t, idx) = dynamic;
		end


		%%Add noise.
		W = randn(N, N);
		if( noiseCorrCoeff > 0 )
			red     = noiseCorrCoeff.*Wprev + sqrt(1-noiseCorrCoeff^2).*W;
			Wprev   = red;
			spatial = spatial + sqrt(noiseVar).*red + noiseMean;
		else	
			spatial = spatial + sqrt(noiseVar).*W + noiseMean;
		end

		Z(:,:,t) = spatial;
	end

	Zmean = mean(Z, 3);
	Z     = Z - repmat(Zmean, [1, 1, T]);
	for i=1:nSources,
		Zmean = mean( squeeze(trueZ(:, :, :, i)), 3);
		trueZ(:, :, :, i) = trueZ(:, :, :, i) - repmat(Zmean, [1, 1, T]);  
	end


	%for t=1:T,
	%	minVal  = prctile(Z(:), 2);
	%	maxVal  = prctile(Z(:), 98);
	%	%minVal  = min(Z(:));
	%	%maxVal  = max(Z(:));
	%	spatial = Z(:, :, t);
	%	imagesc(spatial, [minVal, maxVal]);
	%	colorbar;
	%	pause(.3);
	%
	%end

end
