function plot_PASF_outputs(Z, C, CInfo, SDFInfo, plotPhaseModulus, bfilename)
%
%% AUTHOR    : Mohammad Khabbazian
%% $DATE     : 24-Apr-2017 $
%% $Revision : 1.00 $
%% DEVELOPED : R2016a
%% FILENAME  : pasf.m
%

	EV   = SDFInfo.Eigenvalues; 
	EVec = SDFInfo.Eigenvectors; 
	TR   = SDFInfo.Traces;
	Ph   = SDFInfo.PhaseArray;

	%scrsz = get(groot,'ScreenSize');
	scrsz = [1, 1, 1920, 1080];
	close all;
	set(gcf,'visible','off');
	
	[N, M] = size(EV);
	[d1, d2, d3] = size(Ph);
	

	plot_eigenvalues(EV, TR, bfilename, 1:floor(M/2) );

	plot_clusters(EV, C, CInfo, [d1, d2], bfilename);


	if plotPhaseModulus == false,
		return;
	end

	close all;

	figureHandle = figure('Position', [1 scrsz(4) scrsz(3) scrsz(4)]);
	set(gcf,'visible','off');

	nFrames = floor(size(Ph,3)/2)+1;


	myVideo = VideoWriter( strcat(bfilename,'_EVecsPhaseModulus_VS_Freq', '.avi') );
	myVideo.FrameRate = 4;
	open(myVideo);

	nEVs = size(Ph, 4);
	for i=1:nFrames,
		clf;
		for j=1:nEVs,

			tmp = Ph(:, :, i, j);
			mn  = prctile(tmp(:), 1)  - eps;
			mx  = prctile(tmp(:), 99) + eps;

			subaxis(2, nEVs, j, 'Spacing', 0.01, 'Padding', 0, 'Margin', 0.01, 'SpacingVert', 0.02);
			imagesc( Ph(:, :, i, j), [mn, mx] );
			axis equal tight;
			axis off; colorbar;
			title( strcat('EV', num2str(j), ', EV ', num2str(EV(j,i)), ', C #', num2str(C(j,i)) ) );
			set(gca, 'XTickLabel', '', 'YTickLabel', '')

			subaxis(2, nEVs, nEVs+j, 'Spacing', 0.01, 'Padding', 0, 'Margin', 0.01, 'SpacingVert', 0.02);
			imagesc( abs(EVec(:, :, i, j)) );
			axis equal tight;
			axis off; colorbar;
			title( strcat('F#', num2str(i), ', EV#', num2str(j), ', EV=', num2str(EV(j,i)) ) );
			set(gca, 'XTickLabel', '', 'YTickLabel', '')

		end

		axesHandles = findobj(get(figureHandle,'Children'), 'flat', 'Type', 'axes');
		%axis(axesHandles, 'square')

		
		F = getframe(gcf);
		writeVideo(myVideo, F);
	end
	close(myVideo);

end





function plot_clusters(Eigenvalues, Clusters, ClusterInfo, spatialDim, bfilename)

	[N, M] = size(Eigenvalues);
	freqs  = 1:(floor(M/2)+1);

	figure;
	hold on;
	for i=1:N,
		plot( freqs, Eigenvalues(i, freqs) );
		set ( gca, 'XLim', [min(freqs), max(freqs)] );
	end
	hold off;
	xlabel('Frequency');
	ylabel('Eigenvalue');
	filename = strcat(bfilename, '_topEVsDist');
	set(gcf, 'PaperOrientation', 'landscape');
	print(filename, '-depsc');


	close all;
	subplot(2, 1, 1); 
	hold on;
	for i=1:N,
		plot( freqs, Eigenvalues(i, freqs) );
		set ( gca, 'XLim', [min(freqs), max(freqs)] );
	end
	xlabel('Frequency');
	ylabel('Eigenvalue');

	subplot(2, 1, 2); 
	imagesc( Clusters(:, freqs) );
	box off;
	cMap = [1, 1, 1; parula(max(Clusters(:)))];
	colormap(cMap);
	h = colorbar;
	set(h, 'Position', [.92 0.1 .05 .85])


	filename = strcat(bfilename, '_SDF_properties');
	set(gcf, 'PaperOrientation', 'landscape');
	print(filename, '-dpdf', '-fillpage');

       close all;
       plot_dendrogram(ClusterInfo.Tree, bfilename );

	close all;
	[nTopEVs, NF] = size(Clusters);
	CC    = max(Clusters(:));
	freqs = 1:floor(NF/2);

	%plot_color=['b'; 'r'; 'g'; 'k'; 'y'; 'c'; 'm'];
	plot_color=cool(CC);
	subaxis(nTopEVs+1, 1, 1, 'Spacing', 0.01, 'Padding', 0, 'Margin', 0.03, 'SpacingVert', 0.03);

	hold on;
	[N, M] = size(Eigenvalues);
	for i=1:N,
		plot(freqs, Eigenvalues(i, freqs) );
		set(gca, 'XLim', [min(freqs) max(freqs)] );
		box off;
	end
	hold off;
	xlabel('Frequency');
	ylabel('Eigenvalues');

	for e=1:nTopEVs,
		subaxis(nTopEVs+1, 1, e+1, 'Spacing', 0.01, 'Padding', 0, 'Margin', 0.03, 'SpacingVert', 0.03);

		hold on;

		for c=1:CC, %%NOTE: we ingore zero.
			indices = ( Clusters(e, freqs) == c );
			plot( freqs(indices),  Clusters(e, indices), '.', 'MarkerSize', 18, 'color', plot_color(c,:) );
			set(gca, 'XLim', [0, max(freqs) ], 'YLim', [1, CC]);
			set(gca, 'XTick', [], 'YTick', 1:CC);
			box off;
		end

		hold off;
	end

	filename = strcat(bfilename, '_clusters_eigenvalue_dist.pdf');
	%set(gcf, 'PaperOrientation', 'landscape');
	print(filename, '-dpdf', '-fillpage');

	close all;

	for itr =1:3,
		CC            = max(Clusters(:));
		figureHandle  = figure();
		subplot_ncols = ceil( sqrt(CC) );
		subplot_nrows = ceil( CC / subplot_ncols );

		for i = 1:CC,
			subaxis(subplot_nrows, subplot_ncols, i, 'Spacing', 0, 'Padding', 0.01, 'Margin', 0.04, 'SpacingVert', 0.04);
			if itr==1,
				tmp  = ClusterInfo.Centers(i, :);
			elseif itr == 2,
				tmp  = ClusterInfo.ModCenters(i, :);
			elseif itr == 3,  
				tmp  = ClusterInfo.PhaseCenters(i, :);
			end

			mn = prctile(tmp(:),  1) - eps;
			mx = prctile(tmp(:), 99) + eps;

			imagesc( reshape(tmp, spatialDim(1), spatialDim(2)), [mn, mx] );
			axis equal tight;
			set(gca, 'XTick', [], 'YTick', []);
			colorbar;
			title( strcat('#', num2str(i), ',', 'S:', num2str( ClusterInfo.CSizes(i) ), ... 
				',', 'E:', num2str(ClusterInfo.CEnergy(i)) ) );
		end
		axesHandles = findobj( get(figureHandle, 'Children'), 'flat', 'Type', 'axes');
		%axis(axesHandles,'square')

		if itr == 1, 
			filename = strcat(bfilename, '_cluster_centroid');
		elseif itr == 2,
			filename = strcat(bfilename, '_cluster_modulus');
		elseif itr ==3 
			filename = strcat(bfilename, '_cluster_phase');
		end

		set(gcf, 'PaperOrientation', 'landscape');
		print(filename, '-dpdf', '-fillpage');

		close all;
	end



	close all;
	plot(1:CC, cumsum(ClusterInfo.CEnergy), '--r*', 'MarkerSize', 10);
	set (gca, 'XLim', [0 CC+1], 'YLim', [0 1]);
	filename = strcat(bfilename, '_cluster_cumsum_energy_ratio');
	set(gcf, 'PaperOrientation', 'landscape');
	print(filename, '-dpdf', '-fillpage');

	close all;
	mem = ClusterInfo.Mem;
	seq = 1:length(mem);
	for i = 1:CC,
		subaxis(ceil(CC/2), 2, i, 'Spacing', 0, 'Padding', 0.01, 'Margin', 0.04, 'SpacingVert', 0.04);
		hold on;
		idxs = find(mem==i);
		plot( idxs, ClusterInfo.D(idxs, i), 'k*');
		idxs = find(mem~=i);
		plot( idxs, ClusterInfo.D(idxs, i), 'r*');
		hold off;
		title( strcat( 'C#', num2str(i), ', CS:', num2str(ClusterInfo.CSizes(i)) ) );
		set( gca, 'XTick', [], 'YTick', [0, 1, 2], 'YLim', [0, 2] );
	end

	filename = strcat(bfilename, '_clusters_distance_dist');
	set(gcf, 'PaperOrientation', 'landscape');
	print(filename, '-dpdf', '-fillpage');


	close all;
	for i = 1:CC,
		subaxis(ceil(CC/2), 2, i, 'Spacing', 0, 'Padding', 0.01, 'Margin', 0.04, 'SpacingVert', 0.04);
		hist( ClusterInfo.D(:, i));
		title( strcat( 'C#', num2str(i) ) );
		set( gca, 'XTick', [0, 1, 2]);
	end

	filename = strcat(bfilename, '_clusters_distance_hist');
	set(gcf, 'PaperOrientation', 'landscape');
	print(filename, '-dpdf', '-fillpage');


end



function plot_eigenvalues(EV, Traces, bfilename, frequencies)

	close all;
	[N, M] = size(EV);

	if ~exist( 'frequencies' ),
		frequencies = 1:floor(M/2);
	end

	subplot(ceil( (N+1)/2), 2, 1); 
	plot_EVs( Traces, frequencies, 'Omega', 'Trace' );


	for i=1:N,
		subplot ( ceil( (N+1)/2), 2, i+1 ); 
		plot_EVs( EV(i, :), frequencies, 'Omega', strcat('EV(', num2str(i), ')' ) );
	end

	filename = strcat(bfilename, '_SDF_eigenvalues');
	set(gcf, 'PaperOrientation', 'landscape');
	print(filename, '-dpdf', '-fillpage');

	close all;
	All_EVs = sort(EV(:), 'descend');
	subplot(3, 1, 1);
	plot( 1:length(All_EVs), All_EVs );
	set (gca, 'XLim', [1 length(All_EVs)] );
	%xlabel('Sorted Eigenvalues');

	subplot(3, 1, 2);
	semilogy( 1:length(All_EVs), All_EVs );
	set (gca, 'XLim', [1 length(All_EVs)] );
	%xlabel('Sorted Eigenvalues');

	subplot(3, 1, 3);
	loglog( 1:length(All_EVs), All_EVs );
	set (gca, 'XLim', [1 length(All_EVs)] );
	xlabel('Sorted Eigenvalues');

	filename = strcat(bfilename, '_sorted_eigenvalues');
	set(gcf, 'PaperOrientation', 'landscape');
	print(filename, '-dpdf', '-fillpage');

end




function plot_EVs( ts, freq, xlab, ylab)

	Ymin = -eps; 
	Ymax = max( ts(freq) ) + eps;
	M = length(ts);

	omega = freq./M;
	plot( omega, ts( freq ) );
	set (gca, 'XLim', [min(omega) max(omega)], 'YLim', [Ymin Ymax]);

	xlabel( xlab );
	ylabel( ylab );

end



function plot_dendrogram(Tree, bfilename)

	close all;
	dendrogram(Tree);
	print( strcat(bfilename, '_HC_tree.pdf'), '-dpdf', '-fillpage');

	dendrogram(Tree, size(Tree, 1) );
	print( strcat(bfilename, '_HC_tree_all.pdf'), '-dpdf', '-fillpage');


end


