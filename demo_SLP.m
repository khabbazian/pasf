
load('data/slp.mat');
data = permute(data, [3, 2, 1]);


dirname   = strcat( 'SLP' );
bfilename = strcat( dirname, '/slp' );
mkdir( dirname );

myOpts = struct('cmethod', 'phase', ... 
	'nTopEVs', 2, 'outlierThreshold', 0.7, ...
	'boolParfor', false, ... 
	'saveData', 0, 'boolUseSavedData', 0, ... 
	'kmethod', 'triangle', ...
	'spans', 21, 'errorRate', 0.1, 'bfilename', bfilename);


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
