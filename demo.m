

simOpts=struct('noiseVar', 1, 'noiseCorrCoeff', 0);
data = rotating_energy_sources(500, simOpts);

tic;
Z = pasf(data, 2); 
toc


%%Now, plot the output.
close all;
scrsz = get(groot,'ScreenSize');
figureHandle = figure('Position', [1 scrsz(4) scrsz(3) scrsz(4)]);
set(gca, 'XTickLabel','')

[d1, d2, d3, d4] = size(Z);
minVec = prctile(reshape(Z,d1*d2*d3,d4), 1);
maxVec = prctile(reshape(Z,d1*d2*d3,d4), 99);

for i = 1:d3,
	for j = 1:d4,

		subplot(2, 2, j);
		imagesc( Z(: ,: ,i ,j), [minVec(j), maxVec(j)] );
		colorbar;

		if j == d4,
			title('Input Signal');
		elseif j == d4-1,
			title('Error');
		else
			title( strcat('Signal #', num2str(j)) );
		end

	end

	axesHandles = findobj(get(figureHandle, 'Children'), 'flat', 'Type', 'axes');
	axis(axesHandles,'square')
	pause(0.02);
end

