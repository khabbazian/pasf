
data = propagating_energy_sources();

dirname   = strcat( 'PES' );
bfilename = strcat( dirname, '/pes' );
mkdir( dirname );

myOpts = struct( 'cmethod', 'phase', ...
	'spans', 21, 'detrend', 2, ...
	'errorRate', 0.1, ...
	'boolParfor', false, ... 
	'saveData', 0, 'boolUseSavedData', 0, ...
	'bfilename', bfilename);

tic;
Z = pasf(data, 4, myOpts); 
toc;

[d1, d2, d3, d4] = size(Z);
rng = prctile( Z(:), [1 99] ) + [-eps eps];
for f = 1:d3,
	frame = Z(:,:,f,1);
	for c = 2:d4,
		frame = horzcat(frame, Z(:, :, f, c) );
	end
	imagesc( frame, rng);  
	title('Component #1, Component #2, Component #3, Component #4, Remainder, and Input Signal');
	daspect([1 1 1]);
	axis off;
	pause(0.5);
end
