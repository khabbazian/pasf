
function MVTS=propagating_energy_sources()

	%% Initializing the tuning parameters. 
	N = 20;
	T = 1000; 
	bandwidth = 25; 

	centers = [ [1; 1], [1; N], [N; 1], [N; N] ];
	%centers = [ [1; 1], [1; N] ];
	nSources = size(centers, 2);

	%% This is the variable that the function returns.
	MVTS=zeros(N, N, T); 

	coeffs = [ [0.9; -0.5], [0.9; -0.8], [-0.9; -0.5], [-0.9; -0.8]]; 

	len = 2*N+T+N; 	
	AR = zeros(len, size(centers, 2));
	for idx = 1:nSources,
		noise = randn(1, len);
		AR(:, idx) = filter(1, coeffs(:, idx), noise); % Make AR processes by filtering white noise
	end

	[X,Y] = meshgrid(1:N, 1:N);
	C = [X(:) Y(:)]';

	l1Distances = zeros(N*N, nSources);
	decayRates  = zeros(N*N, nSources);
	for idx=1:nSources,
		differences         = C - repmat(centers(:, idx), 1, N*N);
		decayRates(:, idx)  = exp(-sqrt(sum(differences.^2, 1))/bandwidth);
		l1Distances(:, idx) = sum(abs(differences), 1);
	end

	for t=1:T,
		spatial = zeros(N, N);

		for idx = 1:nSources,
			indeces = (X(:)-1)*N + Y(:); %%TODO double check it.
			spatial(indeces) = spatial(indeces) + AR(t+(2*N-l1Distances(indeces, idx)), idx).*decayRates(indeces, idx);
		end
		MVTS(:, :, t) = spatial;
	end
end



