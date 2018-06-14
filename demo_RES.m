

%% Low level of independent noise
simOpts = struct('noiseVar', 0.16, 'noiseCorrCoeff', 0, 'sourceEnergy', 6.3);

%% Medium level of correlated noise 
%simOpts = struct('noiseVar', 4, 'noiseCorrCoeff', 0.5, 'sourceEnergy', 6.3);
%
%% High level of highly correlated noise
%simOpts = struct('noiseVar', 16, 'noiseCorrCoeff', 0.8, 'sourceEnergy', 6.3);

data    = rotating_energy_sources(70, simOpts);

dirname   = strcat( 'RES-rn' );
bfilename = strcat( dirname, '/res' );
mkdir( dirname );

myOpts = struct( 'cmethod', 'phase', ...
	'boolParfor', false, ... 
	'boolUseSavedData', 0, ...
	'errorRate', 0.1, ...
	'bfilename', bfilename);


tic;
Z = pasf(data, 2, myOpts);
toc;

[d1, d2, d3, d4] = size(Z);
rng = prctile( Z(:), [1 99] ) + [-eps eps];
for f = 1:d3,
	frame = Z(:,:,f,1);
	for c = 2:d4,
		frame = horzcat(frame, Z(:, :, f, c) );
	end
	imagesc( frame, rng);  
	title('Component #1, Component #2, Remainder, and Input Signal');
	daspect([1 1 1]);
	axis off;
	pause(0.5);
end

