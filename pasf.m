function [Z_dynamics, Clusters, Eigenvectors, Eigenvalues, PhaseArray] ... 
		= pasf(Z, nSpatialComponents, options)
% PASF Phased-Aligned Spectral Filtering for Decomposing Spatiotemporal Dynamics.  
% PASF seeks to extract dynamics from spatiotemporal signal via data
% assimilation and modeling.  PASF assumes that the observed spatiotemporal
% data represent superimposed lower-rank smooth oscillations and movements from
% a generative dynamic system, mixed with higher-rank random noises. Separating
% the signals from noises is essential to visualize, model and understand these
% lower-rank dynamic systems. PASF uses a filtering framework for identifying
% lower-rank dynamics and its components that are embedded in a high
% dimensional spatiotemporal system. The approach is based on the structural
% decomposition and phase-aligned construction in the frequency domain.  
% 
% INPUT 'Z' is the spatiotemporal signal of dimension d1 x d2 x T where the
% 		third dimension represents the time samples.
% 	'nSpatialComponents' is the number of dynamic spatial components. PASF
% 		decomposes Z into 'nSpatialComponents' spatial components where
% 		each represents a separate dynamics.
% 	'options' is a struct with optional input parameters. Parameters are:
%		'boolDemeanInput' is a logical variable that indicates if the input
%			should be demeaned before running the PASF. The default is 1.
%		'nTopEVs' is the number of top eigenvectors PASF picks at each
%			frequency from the estimated spectral density function. The
%			default value is nSpatialComponents.
%		'errorRate' is the maximum error rate. The function computes
%			the eigenvalue threshold based on this rate. It should be
%			determine based on SNR.
%		'taperLen' is the taper length that is used for estimating the spectral
%			density function. The default value is zero.
%		'kmethod' is the type of smoothing kernel that is used to estimate
%			the spectral density function. The default is 'triangle'.
%		'spans' is the bandwidth of the smoothing kernel that is used
%			to estimate the spectral density function. The default value is
%			floor(log(size(Z,3))).
%		'detrend' is a variable that indicates if the trend of the signal
%			should be removed for estimating the spectral density
%			function. 0 does not detrend. 1 subtracts the mean. 2
%			subtracts the linear trend. The default is 0.
%		'boolSaveData' is a logical variable that indicates if the
%			PASF's intermediate results should be saved on the hard-drive.
%			It is useful for the large inputs. The default is 0.
%		'bfilename' is the base file name that function uses to save data if
%			indicated. The default is 'pasf_out'.
%		'boolQuietly' is a logical variable that indicates if the progress 
%			method should not be reported. The default is 0. 
% 		
% OUTPUT 'Z_dynamics' is the input spatiotemporal signal that is decomposed
% 	into nSpatialComponents+2 where the last two are noise and the input signal
% 	after demeaning if specified. Z_dynamics is an array of size
% 	d1 x d2 x T x (nSpatialComponents+2). 
%
% Details	
% 	For more details about the PASF refer to https://arxiv.org/abs/1604.04899
% 
% Examples:
% 
% simOpts = struct('noiseVar', 1, 'noiseCorrCoeff', 0);
% data    = rotating_energy_sources(500, simOpts);
% Z       = pasf(data, 2);
% 
% %% Now plot the output 
% set(gca, 'XTickLabel','')
% 
% [d1, d2, d3, d4] = size(Z);
% minVec = prctile(reshape(Z,d1*d2*d3,d4), 1);
% maxVec = prctile(reshape(Z,d1*d2*d3,d4), 99);
% 
% for i = 1:d3,
% 	for j = 1:d4,
% 		subplot(2, 2, j);
% 		imagesc( Z(: ,: ,i ,j), [minVec(j), maxVec(j)] );
% 		colorbar;
% 
% 		if j == d4,
% 			title('Input Signal');
% 		elseif j == d4-1,
% 			title('Error');
% 		else
% 			title( strcat('Signal #', num2str(j)) );
% 		end
% 	end
% 
% 	axesHandles = findobj(get(figureHandle, 'Children'), 'flat', 'Type', 'axes');
% 	axis(axesHandles,'square')
% 	pause(0.02);
% end % 

 
if nargin == 1
	nSpatialComponents = 1;
end

%% setting up the default values
pasfOpts = struct('nTopEVs', nSpatialComponents, ...
'errorRate', 0.01, ...
'taperLen', 0, 'kmethod', 'triangle', 'spans', floor(log(size(Z,3))), ...
'boolDemeanInput', 1, 'boolQuietly', 0, ...
'detrend', 0, 'boolZeroPad', 0, 'boolSaveData', 0, 'bfilename', 'pasf_out');


if nargin > 2
	fields = fieldnames(options);
	for f=1:length(fields),
		pasfOpts.(fields{f}) = options.(fields{f});
	end
end


if pasfOpts.boolDemeanInput == 1,
	Trend = mean(Z, [3]);
	Z = Z - repmat(Trend, [1, 1, size(Z,3)]);
end


[d1, d2, d3] = size(Z);
if(d3 > 1)
	Z = permute(reshape(Z, d1*d2, d3), [2,1]); %NOTE: each column is a time series.
end

nTimes   = size(Z, 1);
nSpatial = size(Z, 2);
nTopEVs  = pasfOpts.nTopEVs;  

if ~pasfOpts.boolQuietly,
	disp(pasfOpts);
	disp('done with the initial steps.');
end


[PVec_FT_Cxx, PPhase_FT_Cxx, Eigenvalues, NF, traceSums] = mvspec(Z, pasfOpts);

if ~pasfOpts.boolQuietly,
	disp('done with estimating the spectral density matrices and their eigen decomposition.');
end


if pasfOpts.boolSaveData,
	filename = strcat(bfilename, '_intermediate_data.mat');
	save(filename);
end
%load(filename);

%%NOTE: uncomment the following to compile the c code.
if exist('unwrap2D') ~= 3,
	mex unwrap2D.c
end

%scoreMask = Eigenvalues > evThreshold;

sortedEVs   = sort(Eigenvalues(:), 'descend');
indicesMask = cumsum(sortedEVs)./traceSums > 1 - pasfOpts.errorRate;
if(sum(indicesMask)>0)
	eigenValueThreshold = max( sortedEVs(indicesMask) );
else
	eigenValueThreshold = 0;
	warning('set the Delta(the eigenvalue threshold) to zero,... you may want to increase number of top eigenvalues.');
end

if ~pasfOpts.boolQuietly,
	disp(strcat('eigenvalue threshold is: ', num2str(eigenValueThreshold), '.') );
end

scoreMask = Eigenvalues > eigenValueThreshold;

PPhase_FT_Cxx = reshape(PPhase_FT_Cxx, nSpatial, nTopEVs * NF )';
for i=1:size(PPhase_FT_Cxx, 1)
	PPhase_FT_Cxx(i, :) = reshape( unwrap2D(reshape(PPhase_FT_Cxx(i, :), d1, d2), zeros(d1, d2)), 1, d1*d2);
	%NOTE: subtracting the mean doesn't change the result but makes the plots prettier. 
	PPhase_FT_Cxx(i, :) = PPhase_FT_Cxx(i, :) - mean(PPhase_FT_Cxx(i, :)); 	
end

%%NOTE: keep the unwrapped phase for ploting.
PhaseArray = zeros(d1, d2, NF, nTopEVs);
for j=1:NF, 
	for k=1:nTopEVs,
		PhaseArray(:, :, j, k) = reshape(PPhase_FT_Cxx( (j-1)*nTopEVs+k, :), d1, d2);
	end
end

if ~pasfOpts.boolQuietly,
	disp('done with the phase unwrapping.')
end


%Clusts = zeros(size(PPhase_FT_Cxx, 1), 1); 
Clusts = zeros(size(scoreMask)); 

%%NOTE: kmeans/linkage accepts the point cloud as a matrix whose rows are the data point (NOT columns).
%kmeans with correlation as the similarity metric.
half = (2+ceil((NF-1)/2));
scoreMask(:, half:end) = 0;
scoreMask = scoreMask(:);
[memberships, centers, sumD] = kmeans(PPhase_FT_Cxx(scoreMask(:), :), nSpatialComponents, 'Distance', 'correlation', 'Replicates', 128);
Clusts(scoreMask) = round( memberships );
Clusts(:, half:end) = Clusts(:, half-2:-1:2);
Clusts = Clusts(:); 

%Agglomerative hierarchical clustering with the Ward's method and correlation as the similarity metric.
%UTree = linkage(PPhase_FT_Cxx(scoreMask, :), 'ward', 'correlation');
%Clusts(scoreMask) = cluster(UTree, 'maxclust', nSpatialComponents);

if ~pasfOpts.boolQuietly,
	disp('done with the clustering.');
end

C_omega = PVec_FT_Cxx;
%TODO: Do we need to store B_omega explicitly?
B_omega = zeros(nTopEVs, nSpatial, size(PVec_FT_Cxx, 3));
for i=1:size(PVec_FT_Cxx, 3),
	B_omega(:, :, i) = PVec_FT_Cxx(:, :, i)';
end

Z_dynamics = zeros(nTimes, nSpatial, nSpatialComponents);
for cNum = 1:nSpatialComponents,
	boolClust = (Clusts == cNum);
	mask = repmat( reshape(boolClust , nTopEVs, 1, size(PVec_FT_Cxx, 3)), 1, nSpatial, 1);

	X = zeros(NF, nTopEVs);
	tZ = zeros(size(Z, 1), 1);

	for i=1:size(B_omega,  1),
		for j=1:size(B_omega, 2),
			tZ(1:nTimes) = squeeze(Z(:,j));
			X(:, i) = X(:, i) + squeeze(B_omega(i, j, :).*mask(i, j, :)) .* fft(tZ);
		end
	end

	mask  = permute(mask,[2,1,3]);
	Z_hat = zeros(NF, nSpatial);

	for i=1:size(C_omega, 1),
		for j=1:size(C_omega, 2),
			Z_hat(:, i) = Z_hat(:, i) + ifft(squeeze(C_omega(i, j, :).*mask(i, j, :)) .* squeeze(X(:, j)) );
		end
	end
	Z_dynamics(:, :, cNum) = real(Z_hat(1:nTimes, :)); 
end

if ~pasfOpts.boolQuietly,
	disp('done with the filtering/separating dynamic components.');
end


%%preparing the outputs
Z_dynamics = permute(Z_dynamics, [2,1,4,3]);
Z_dynamics = reshape(Z_dynamics, d1, d2, d3, nSpatialComponents);
Z          = reshape(permute(Z, [2, 1]), d1, d2, d3);
Err        = Z - sum(Z_dynamics, 4);
Z_dynamics(:, :, :, nSpatialComponents + 1) = Err; 
Z_dynamics(:, :, :, nSpatialComponents + 2) = Z;

Clusters     = reshape(Clusts, nTopEVs, NF);
%NOTE: PVec_FT_Cxx \in C^{S x nTopEVs x T};
Eigenvectors = permute(reshape(abs(PVec_FT_Cxx), d1, d2, nTopEVs, NF), [1, 2, 4, 3]);

if ~pasfOpts.boolQuietly,
	disp( strcat('Error rate is: ', num2str(sum(Err(:).^2)/sum(Z(:).^2)), '.'));
end

end %end of pasf


%%It estimates the multivariate spectral function.
function [PVec_FT_Cxx, PPhase_FT_Cxx, Eigenvalues, N, traceSums] = ...
		mvspec(X, opts)

N0 = size(X,1); 
N  = N0;
ncols  = size(X,2);
xfreq  = 1;
nTopEVs = opts.nTopEVs;

if(opts.detrend == 1 )
	X = detrend(X);
else 
	if( opts.detrend == 2 )
		X = detrend(X, 'constant');
	end
end

if opts.taperLen > 0,
	w = taper(N, opt.taperLen)';
	X = X .* repmat(w, 1, ncols);
end


if opts.boolZeroPad ==  1, 
	X     = [X; zeros(N-1, ncols)];
	N     = size(X,1); 
	ncols = size(X,2);
end

xfft = fft(X);
clearvars X; %we don't need it anymore.

if ~opts.boolQuietly,
	disp('done with xfft computation');
end


spans = 2*floor(opts.spans/2)+1; %%just to make the smoothing easier
switch opts.kmethod
	case 'box'
		g = ones(spans,1);
	case 'gaussian'
		g = gausswin(spans);
	case 'triangle'
		g = triang(spans);
	otherwise
		error('kmethod is undefined [box|gaussian|triangle]');
end
g = g/sum(g);

PVec_FT_Cxx   = zeros(ncols, nTopEVs, N); 
PPhase_FT_Cxx = zeros(ncols, nTopEVs, N); 
%PScores_Cxx   = zeros(nTopEVs, N);
Eigenvalues   = zeros(nTopEVs, N);
empericalThreshold = zeros(N,1);

traceSums = 0;
half = ceil(spans/2);
for i=2:floor(N/2)+1, 
%for i=2:N, 
	M = zeros(ncols, ncols);
	%% for 128x128 the following takes 14 sec 
	%NOTE: estimate the spectral density matirx by smoothing the FFT of empirical autocovariance.
	%NOTE: since all the coefficients are positive M is positive semi-definite.
	for j=1:spans, %with the assumption that span is an even number
		M = M + g(j) * get_periodogram(xfft, N, ncols, N0*xfreq, i, j-half);
	end
	traceSums = traceSums + 2*trace(M);
	%M = 1/N .* M;
	%% for 128x128 with nTopEVs=4 the following takes 12 sec 
	[V, D] = svds(M, nTopEVs, 'L'); 
	D = diag(D);

	PVec_FT_Cxx(:, :, i) = V; 
	PVec_FT_Cxx(:, :, N-i+2) = conj(V); 
	%NOTE: all the eigenvalues must be positive, so trace(M)==\sum_i \lambda_i. TODO: fix the following accordingly.
	%const = (M(:)'*M(:));
	%PScores_Cxx(:, i)    = (D.^2)/const; %denominator is eqv to trace(M*M):but faster to compute.
	Eigenvalues(:, i)    = D; %NOTE: all the eigenvalues must be positive so abs is unnecessary. 
	Eigenvalues(:, N-i+2)    = D; %NOTE: all the eigenvalues must be positive so abs is unnecessary. 

	%empericalThreshold(i) = median( empericalEValueThreshold(M, 10) );
	PPhase_FT_Cxx(:, :, i) = angle(V);
end

end 


%TODO: change the name of the following function
function pgram = get_periodogram(xfft, N, ncols, const, i, offset)
	idx = i + offset;
	if( idx < 2 ) %ignoring zero freq;
		idx = N + idx - 1; %-1 to ignore zero freq 
	end
	if( idx > N )
		idx = idx - N;
	end
	pgram = zeros(ncols, ncols);
	%pgram = xfft(idx, :).' * conj( xfft(idx, :) )/(const);
	pgram = kron(xfft(idx, :).', conj(xfft(idx, :)))/(const); % eqv to the above line but faster!
end


function w = taper(n, p) 
	if( p > 0.5 || p < 0 )
		error('p must be between 0 and 0.5');
	end

	m = floor(n * p);
	w = 0.5 * (1 - cos(pi * (1:2:2*m-1)/(2*m)) );
	w = [w, ones(1, n - 2*m), flip(w)];
end


%NOTE: This function tries to form a random matrix out of M entries 
% to find a threshold for eigenvalues of the original matrix with a structure.
function topEValue = empericalEValueThreshold(M, nIterations)

	N = size(M,1);

	%TODO subtract the mean of each row from the M and after resampling 
	%add it back. Check the lit to see if someone has developed and studied the bootstrap procedure
	topEValue = zeros(nIterations, 1);
	for i=1:nIterations,
		Mp = reshape(datasample(M(:), N*N), N, N);
		Mp = triu(Mp) + triu(Mp, +1)';
		topEValue(i) = abs(eigs(Mp, 1));
	end

end
