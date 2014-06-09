function varargout = vg_grabCut3D( C, T, mode, varargin )

% vg_grabCut3D : Computes a volume segmentation (vg in the beginning stands for Varun Gulshan ;)
%
% C -> m x n x c x t image volume: m x n is size of each frame, c = number of channels t = number of frames
%                    currently c should be 3 for this code to work
% T -> m x n x t input trimaps.'uint8'
%
% mode: currently following modes are supported (and the corresponding varargin are described
%
%  'single' [default mode] nema_grabCut3D(C,T,'single',opts) 
%   where opts is the structure of options (we have defaults for all options, so partial or even no structure is fine)
%   See later for what are the various elements of opts
%   
%   output is: [s] s= m x n x t 3D segmentation mask (0=background, 255=foreground)
%
%  T is defined as follows (similar to nema_matte_grabcut):
%
%      T is a trimap with the following pixel values:
%         0   - definite background , clamped to be always background
%         64  - definite background used for colour distribution, clamped to be always background
%         96  - background used for color model initialization, not clamped
%         128 - unknown region
%         160 - foreground used for color model initialization, not clamped
%         192 - definite foreground used for colour distribution, clamped to be always foreground
%         255 - definite foreground, clamped to be always foreground
% 
%   opts has the following fields
%
%   gammaTime (50): weighing the binary term across frames with the unary term
%   fgModelTimeWindow(1): the number of frames across which the fg model for a particular frame is learnt.
%                      Has to be odd, because it is symettric about the frame. For the boundary frames,
%                      it will shift the window in the other direction in case this window overflows the boundary
%                      If you give even, it will be rounded to the next odd number
%   bgModelTimeWindow(1): similar to fgModelTimeWindow
%  
%   And also the fields below (these are the same as in nema_matte_grabCut)
%       nmixtures (5): maximum number of mixture components for
%         modelling the foreground and background colours.
%       N (40): width of initial background region.
%       gamma (50): region weighting.
%       bwmorph (1): apply morphological operations after convergence
%         to clean the alpha matte.
%       lambda1 (50): energy coefficient 1.
%       lambda2 (1e3): energy coefficient 2.

%       colour_model ('gmm_bs'): colour model estimator to use.  Possible
%         values include:
%           'floNoWindow' - our texture based model, need to supply additional params for this, described below
%           'gmm_bs' - binary split algorithm.
%           'gmm_var' - variational modelling.
%           'lut' - lookup table.  Note that the lookup table
%             violates the energy minimization scheme in that there
%             is no guarantee that the energy will be minimized at
%             each iteration.
%       gmmf ([]): initial foreground Gaussian mixture model.
%         Doesn't use the trimap to estimate the model.
%       gmmb ([]): initial background Gaussian mixture model.
%         Doesn't use the trimap to estimate the model.
%       debug (0): output debug information.  Set to -1 for quite mode.
%
%    Texture energy options(in case color model is floNoWindow):
%    Provide them in the form opts.texture.parms.<paramName>
%    if color_model is floNoWindow, info in texture.parms needs to be
%        parms.patchSize, parms.preScale, parms.normType (0=no normalisation, 1=weber normalisation of all 3 channels)
%        parms.ss_cluster, parms.numTextons, parms.clusterIter,
%        parms.unary='chiSq'|'logProb' (this is form of unary energy used while calculating it at a particular pixel)
%

if(~exist('mode','var'))
 mode='single';
end;

switch(mode)
case {'single'}
  opts = vgg_argparse( struct( 'nmixtures', 5, 'gamma', 50, ...
            'gammaTime',50, 'fgModelTimeWindow', 1, 'bgModelTimeWindow',1, ...
			      'bwmorph', 1, 'border', 7, ...
			      'lambda1', 50, 'lambda2', 1e3, ...
			      'debug', 0, ...
			      'colour_model', 'gmm_bs', ...
			      'lutBlurSigma', 9/255), ...
		        varargin );
  otherwise,
    error('Unsupported mode %s passed to function\n',mode);
end;

[ nrows, ncols, ndims, nframes ] = size( C );
C = im2double( C );

if(mod(opts.fgModelTimeWindow,2)==0 | mod(opts.bgModelTimeWindow,2)==0),
  error('the time window needs to be odd\n');
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Segment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch(mode)
  case 'single'
    [A]=vg_segment3D(C,T,mode,opts);
end;

varargout{1}=A;
return;

function varargout = vg_segment3D( C, T, mode, varargin )
%% Segment the image volume.
%
% Input parameters:
%   C - Input image.
%   T - Trimap image.
%   mode - mode in which grab cut is running
%   varargin - parameterd depending on the mode
%
% Output parameters(depending on the mode):

switch(mode)
 case {'single'}
   opts=varargin{1};
 otherwise
   error('Unknown mode passed\n');
end;

% ----- Preprocessing (independent of trimap) -----
switch(mode)
  case{'single'}

    switch(opts.colour_model)
      case 'floNoWindow',
        opts.texture.parms.ss_x=1;
        opts.texture.parms.ss_y=1;
        pp.txtLabels=learnLabelsTexture3D(C,opts.texture.parms);
      case 'gmm_bs',
        % Nothing to preprocess
        fprintf('');
      otherwise,
        error('Color model = %s not supported right now\n',opts.colour_model);
    end;

    % --- set up the graph's for minimization --
    if opts.debug > 0
      disp( 'Computing the graph ...' );
    end

    [ pp.E_n, pp.E_w, pp.K_threshold ] = vg_setup_graph3D( C, opts.gamma, opts.gammaTime );
end;


% ---- Pre processing done, all preprocessed info in the structure pp ---


% --------- Performing first iteration ---------------

% If the user has not specified any fg or bg pixels then give error
if(isempty(find(T==192|T==160))), 
  error('No foreground pixels specified, cant proceed\n');
elseif (isempty(find(T==64|T==96))), 
  error('No background pixels specified, cant proceed\n');
end;


% --- Create initial colour models

% Grabcut splits the data matrix into opts.nmixtures clusters by
% dividing the data equally (spatially) along its principal vector.
% We are going to do something similar, but slightly more
% sophisticated.  We use the binary split algorithm for colour
% indexing.

if opts.debug > 0
  disp( 'Estimating initial colour models ...' );
end

switch opts.colour_model
 case 'gmm_bs' % calculate the GMMs using a binary split algorithm
  state.gmmf = ...
    namg_init_gmm3D( C , T==160|T==192, opts.nmixtures, opts.fgModelTimeWindow );
  state.gmmb = ...
    namg_init_gmm3D( C , T==64|T==96, opts.nmixtures, opts.bgModelTimeWindow);
 case 'floNoWindow'
   fgMask=(T==160|T==192);
   bgMask=(T==64|T==96);
   [fgModel,bgModel]=learnFGBGTextures3D(opts.texture, pp.txtLabels, fgMask, bgMask, fgModelTimeWindow, bgModelTimeWindow);
   state.txtFgModel=fgModel;
   state.txtBgModel=bgModel;
   clear fgModel bgModel;
  otherwise
  error( sprintf( 'No such colour model "%s" supported yet', opts.colour_model ) );
end

% --- Now T==160 and T==96 have been used to learn models, convert them into unknown(128)
T(find(T==160|T==96))=128;

if opts.debug > 0
  disp( 'Solving ...' );
end

energy = [];
flow = [];

% --- INITIALIZE D Graph Cut -----
[nrows,ncols,ndims,nframes]=size(C);
state.dgcHandle=mexDGC('initialize',uint32(nrows*ncols*nframes),uint32(pp.E_n));
nEnergies=zeros(4,size(pp.E_w,2));
nEnergies([2 3],:)=pp.E_w;
global scale;
scale=50;

% The command below can be improved by actually adding an editweights command to mexDGC, currently not implemented
[x,x]=mexDGC('minimizeDynamic',state.dgcHandle,uint32([]),int32([]),uint32(pp.E_n),int32(round(scale*nEnergies)));


% --- Now get first segmentation using these models ----
switch opts.colour_model
  case 'floNoWindow'
    [STE]=getFloNoWindowSTE_3D(C,T,state.fgModel,state.bgModel,pp.K_threshold);
    [A,energy(end+1),flow(end+1)]= vg_dgraphcut(state.dgcHandle, STE, nrows,ncols,nframes);

  case {'gmm_bs', 'gmm_var'}
    
    [STE,kf,kb]=getGmm_STE3D(C,T,state.gmmf,state.gmmb,pp.K_threshold);
    [A,energy(end+1),flow(end+1)]= vg_dgraphcut(state.dgcHandle, STE, nrows,ncols,nframes);

end

if(opts.debug>1)
  showEnergies3D(opts,state,C,STE,size(T));
end;


% ---  iterative minimisation
iteration = 1;
dA = 1;
while dA & iteration <= 10
  A0 = A; % save old solution

  switch opts.colour_model
   case {'gmm_bs', 'gmm_var'}
    % assign GMM components to pixels

    %keyboard;
    [state.gmmf,state.gmmb] = updateGmm3D(C,A,T,opts,kf,kb,state.gmmf,state.gmmb);
    %keyboard;
    [STE,kf,kb]=getGmm_STE3D(C,T,state.gmmf,state.gmmb,pp.K_threshold);
    %% estimate segmentation using min cut
    [A,energy(end+1),flow(end+1)]= vg_dgraphcut(state.dgcHandle, STE, nrows,ncols,nframes);
    %flow
    %keyboard;

    case 'floNoWindow',
      fgMask= (A==1 & ( T == 128 | T == 192));
      bgMask= (A==0 & ( T == 128 | T == 64 ));
      [fgModel,bgModel]=learnFGBGTextures3D(opts.texture, pp.txtLabels, fgMask, bgMask, fgModelTimeWindow, bgModelTimeWindow);
      state.txtFgModel=fgModel;
      state.txtBgModel=bgModel;
      clear fgModel bgModel;
    
      STE=getFloNoWindowSTE_3D(C,T,state.fgModel,state.bgModel,pp.K_threshold);
      [A,energy(end+1),flow(end+1)]= vg_dgraphcut(state.dgcHandle, STE, nrows,ncols,nframes);

  end
      
  %flow % Temporary, just to check if everything is correct
  if iteration & (flow( end ) > flow( end-1 ))
    warning( 'Increase in energy detected' );
  end
  
  dA = any( A0(:) ~= A(:) );
  if opts.debug > 0
    disp( sprintf( 'Iteration: %d; Energy: %e', iteration, energy( end ) ) );
  end
  
  iteration = iteration + 1;
  
  if(opts.debug>1)
    showEnergies3D(opts,state,C,STE,size(T));
  end;

end % while dA

switch(mode)
  case 'single'
    varargout{1}=A;
    mexDGC('cleanup',state.dgcHandle);
end;

function [E_n, E_w, K] = vg_setup_graph3D(C,gamma,gammaTime)
%% Compute the edge costs.
%
% Input parameters:
%   C - Input image volume
%   gamma - wieghing of binary terms in one frame
%   gammaTime - weighing of binary terms across frames
%
% Output parameters:
%   E_n - Pairwise node labels.
%   E_w - Pairwise terms.
%   K - Large number used to fix nodes in unary energies.
%
%   Note: The nodes are given numbers in a column majoring order for each frame
%   So if there are n pixels in each frame, 1st frame nodes are numbers 1 to n, 2nd frame is n+1 to 2n and so on...

beta = namg_beta3D( C ); % A seperate beta is computed for each frame
betaTime = namg_betaTime3D( C ); % A seperate beta is computed for each frame pair

[ nrows, ncols, ndims,nframes ] = size( C );
N = nrows * ncols; % the number of nodes in one frame

% 8-connected neighbourhood
roffset = [ 1, 1, 1, 0 ];
coffset = [ -1, 0, 1, 1 ];
nbrhoodSizeFrame=length(roffset);
nbrhoodSizeTime=2;
% Don't specify edges twice!
%roffset = [ 1, 1, 1, 0, -1, -1, -1, 0 ];
%coffset = [ -1, 0, 1, 1, 1, 0, -1, -1 ];

%d = sqrt( roffset.^2 + coffset.^2 );

r = 1:nrows;
c = 1:ncols;
[r, c] = meshgrid( r, c );

r = repmat( r(:)', [1, length( roffset )] );
c = repmat( c(:)', [1, length( coffset )] );
roffset = repmat( roffset, [N, 1] );
coffset = repmat( coffset, [N, 1] );
roffset = roffset(:)';
coffset = coffset(:)';
rs = r + roffset;
cs = c + coffset;
ds = sqrt( roffset.^2 + coffset.^2 );
dsTime = ones(1,N);

pixel = r + (c-1)*nrows;
pixels = rs + (cs-1)*nrows;

valid = find( rs >= 1 & rs <= nrows & cs >= 1 & cs <= ncols );
E_nFrameTemplate = [pixel( valid );
                    pixels( valid )];
E_nTimeTemplate = [ (1:N) ; (1:N)+N];
edgesPerFrame=size(E_nFrameTemplate,2);
E_n=zeros(2,nframes*edgesPerFrame+(nframes-1)*N);
E_w=zeros(2,nframes*edgesPerFrame+(nframes-1)*N);

% Set up the edges in each frame first
for i=1:nframes
  curFrame=C(:,:,:,i);
  Cc = nema_impixel( curFrame, c(valid), r(valid) )';
  Cb = nema_impixel( curFrame, cs(valid), rs(valid) )';
  E_n(:,1+(i-1)*edgesPerFrame:i*edgesPerFrame)=(i-1)*N+E_nFrameTemplate;
  E_w(:,1+(i-1)*edgesPerFrame:i*edgesPerFrame) = repmat( gamma * exp( - beta(i) * sum( ( Cc - Cb ).^2, 1 ) ) ./ ...
	        ds(valid), [ 2, 1 ] );
end

% Set up the edges across frames now
for i=1:(nframes-1)
  curFrame=C(:,:,:,i);
  nxtFrame=C(:,:,:,i+1);
  Cc = nema_impixel( curFrame)';
  Cb = nema_impixel( nxtFrame)';
  stIndex=nframes*edgesPerFrame+(i-1)*N+1;
  endIndex=stIndex+N-1;;
  E_n(:,stIndex:endIndex)=(i-1)*N+E_nTimeTemplate;
  E_w(:,stIndex:endIndex) = repmat( gammaTime * exp( - betaTime(i) * sum( ( Cc - Cb ).^2, 1 ) ) ./ ...
	        dsTime, [ 2, 1 ] );
end;

% Modifying nick's way of finding K
% Nick's way was to take the K=maximum over each pixel of (sum of pairwise energies of that pixel)
% Nick finds out this sum of pairwise energy for each pixel, and finds the exact K
% we can just ease some calculation and use an upper bound = max size of clique * maximum pairwise energy

K=(nbrhoodSizeFrame+nbrhoodSizeTime)*max(E_w(1,:));

% Might be wise to check this calculation of K
%blah = zeros( N, length( roffset ) / N );
%blah( valid ) = E_w( 1, : );
%blah = sum( blah, 2 );
%K = 1 + max( blah );

function gmm = namg_init_gmm3D( C, mask, nmixtures, timeWindow )
%% Initialize a GMM for each frame via vector quantization.
% 
% Input parameters:
%   C - Input image volume (m x n x 3 x t)
%   mask  - m x n x t (a 3D mask indication which pixels to consider while building the gmm's
%   nmixtures - The number of mixtures.
%   timeWindow - the frames to be considered around while building model for the current frame
%
% Output parameters:
%   gmm - A Gaussian mixture model structure.
[rows,cols,channels,frames]=size(C);
if(frames < timeWindow)
  error('Time window = %d is greater than the number of frames = %d\n',timeWindow,frames);
end;

halfTimeWindow=floor(timeWindow/2);

indexBeginFrame=zeros(frames+1);
indexBeginFrame(1)=1;
for i=1:frames,
  inds{i}=find(mask(:,:,i)); % inds{i} has the pixel indices for the ith frame, from which to build the gmm model
  indexBeginFrame(i+1)=indexBeginFrame(i)+length(inds{i});
end;

Cnew=zeros(channels,indexBeginFrame(frames+1)-1);
for i=1:frames,
  Cnew(:,indexBeginFrame(i):indexBeginFrame(i+1)-1)=nema_impixel(C(:,:,:,i),inds{i})';
end;

% Compute gmm models for frames where the window fits into the frames
for j=(1+halfTimeWindow):(frames-halfTimeWindow),
  jStart=(j-halfTimeWindow);
  jEnd=(j+halfTimeWindow);
  Ctmp=Cnew(:,indexBeginFrame(jStart):indexBeginFrame(jEnd+1)-1);
  if(size(Ctmp,2))==0, error('Cant initialize GMM model for frame number %d, No input colors\n',j); end;

  [cluster, colours] = nema_vector_quantize( Ctmp, nmixtures );
  for i = 1:nmixtures
    idx = find( cluster == i );
    gmm(j).mu( :, i ) = mean( Ctmp( :, idx ), 2 );
    if(length(idx)<=2),
      gmm(j).sigma(:,:,i) = 0.00001*eye(3);
    else
      gmm(j).sigma( :, :, i ) = cov( Ctmp( :, idx )' );
    end
    gmm(j).pi( i ) = length( idx );
  end
  gmm(j).pi = gmm(j).pi / sum( gmm(j).pi );
end

% Copy the gmm models from the middle frames to the boundary frames, as for the boundary frames
% the window went outside the frames provided
for j=1:halfTimeWindow, 
  gmm(j)=gmm(1+halfTimeWindow);
end;
for j=(frames-halfTimeWindow+1):frames
  gmm(j)=gmm(frames-halfTimeWindow);
end;

function [STE,kf,kb]=getGmm_STE3D(C,T,gmmf,gmmb,K)
% Find the Unary energies , given gaussians for each frame
% Also return which gaussian each pixel in each frame was assigned to, used when updating the gaussians
% Input parameters:
%   C - Input image volume
%   gmm - Gaussian mixture models for each frame
%   T - Trimap image.
%   K - the unary energy to assign when a pixel is hard labeled at fg/bg
%
% Output paramters:
%   STE - the unary energies in one big 2xN array (N=number of frames*number of pixels in each frame)

[ nrows, ncols, ndims, nframes ] = size( C );

nPerFrame=nrows*ncols;

kf = zeros( nrows, ncols, nframes, 'uint8' );
kb = zeros( nrows, ncols, nframes, 'uint8' );
STE= zeros(2,nrows*ncols*nframes);

for i=1:nframes,
  computeIndices=find(T(:,:,i)==128|T(:,:,i)==64|T(:,:,i)==192);
  [steB,kb(:,:,i)]=vg_assignGmm_computeE(C(:,:,:,i),computeIndices,gmmb(i));
  [steF,kf(:,:,i)]=vg_assignGmm_computeE(C(:,:,:,i),computeIndices,gmmf(i));

  idx = find( T(:,:,i) > 128 );
  steB(idx ) = K;
  steF(idx ) = 0; 
  idx = find( T(:,:,i) < 128 );
  steB(idx ) = 0;
  steF(idx ) = K;

  STE(1,1+(i-1)*nPerFrame:i*nPerFrame)=steB;
  STE(2,1+(i-1)*nPerFrame:i*nPerFrame)=steF;

end

function [E,k]=vg_assignGmm_computeE(C,indices,gmm)
% Computes which gmm a pixel is assigned to and its corresponding negative log probabilily also
% Inputs:
%   C -> m x n x dim image
%   indices -> indcies where to compute these gmm assigments and corresponding prob. You dont want to do it at all
%              pixels because some pixels are just useless, they are hard labeled, and no color models are learnt from them
%   gmm -> the gmm model
%
% Outputs:
% E -> the negative log probabilities for each pixel specified in indices, (0 for others)
% k -> the corresponding gaussian which maximises the probability

[ nrows, ncols, ndims ] = size( C );
if length( gmm.pi ) < 256
  k = zeros( nrows, ncols, 'uint8' );
  E = zeros(1,nrows*ncols);
else
  error( 'The number of clusters exceeds 256' );
end

%indices = find( T == 128 | T == 64 | T == 192 );
pixels = nema_impixel( C, indices )';
npixels = length( indices ); %nrows * ncols;
logp = zeros( length( gmm.pi ), npixels );

for i = 1:length( gmm.pi )
  logp( i, : ) = log( gmm.pi( i ) ) + ...
      nema_lognorm( pixels, gmm.mu( :, i ), gmm.sigma( :, :, i ) );
end

[maxp, maxpi] = max( logp, [], 1 );
k( indices ) = maxpi;
E( indices) = -maxp;  % It negative log probability, so putting the - sign 

function [gmmf,gmmb] = updateGmm3D(C,A,T,opts,kf,kb,curGmmf,curGmmb)
% Fit foreground and background gmm's, given the segmentation A
%
% Input parameters:
%   C - Input image volume
%   A - Binary segmentation volume
%   T - trimap
%   opts - the entire opts structure, though we use only a few fields here, but just to keep the # of arguments small
%   kf - Gaussian labels of pixels for each frame for its corr fg gaussian
%   kb - Gaussian labels of pixels for each frame for its corr bg gaussian
%   curGmmf,curGmmfb - the current gmm models. required because you need the find the gaussian labels for the window
%   and not just the current frame

[ nrows, ncols, ndims, nframes ] = size( C );

% Doing for FG first
halfTimeWindow=floor(opts.fgModelTimeWindow/2);
kfWindow=zeros(nrows,ncols,opts.fgModelTimeWindow,'uint8');

% i = corr frame for which we are fitting the new gmm
for i=(1+halfTimeWindow):(nframes-halfTimeWindow),
  iStart=i-halfTimeWindow;
  iEnd=i+halfTimeWindow;
  % iStart and iEnd are the range of frames to which we fit our gmm's
  % kfWindow is the labels given by the current fg gaussian to the pixels in this window
  for j=iStart:iEnd,
    % j indexes the frames of this window, to fill in kfWindow
    if(isequal(curGmmf(j),curGmmf(i))),
      kfWindow(:,:,1+j-iStart)=kf(:,:,j);
    else,
      %computeIndices=find(T(:,:,j)=128|T(:,:,j)=64|T(:,:,j)=192);
      computeIndices=find(T(:,:,j)==128|T(:,:,j)==192); % No need to consider T=64, as it is hard background
      [xNotReq,kfWindow(:,:,1+j-iStart)]=vg_assignGmm_computeE(C(:,:,:,j),computeIndices,curGmmf(i));
      clear xNotReq;
    end
  end

  % So now we have kfWindow containing all the frames in the window with foreground gausssian labels

  for iMix=1:opts.nmixtures
    CiMix=[];
    for j=1:opts.fgModelTimeWindow
      iFrame=i-halfTimeWindow+(j-1);
      idxi = find( kfWindow(:,:,j) == iMix & A(:,:,iFrame)==1 & (T(:,:,iFrame)==128|T(:,:,iFrame)==192)  );
      CiMix = [CiMix; nema_impixel( C(:,:,:,iFrame), idxi )];
    end
    if (size(CiMix,1)~=0),
      gmmf(i).mu( :, iMix ) = mean( CiMix, 1 )';
      if(size(CiMix,1)<=2)
        gmmf(i).sigma(:,:,iMix)=0.00001*eye(3);
      else
        gmmf(i).sigma( :, :, iMix ) = cov( CiMix );% + eye( 3 ) * 25/255/255; % add noise
      end
      % If the cov-matrix is badly scaled, add some noise. It can be badly scaled for 2 reasons:
      % <= 2 members got assigned to this gauusian
      % the covariance is actually very low, because the pixels correspond to a uniform region
      if(det(gmmf(i).sigma(:,:,iMix))<1e-20), gmmf(i).sigma(:,:,iMix)=gmmf(i).sigma(:,:,iMix)+0.00001*eye(3); end;
      gmmf(i).pi( iMix ) = size(CiMix,1);
    else
      gmmf(i)=curGmmf(i);
      gmmf(i).pi( iMix ) = 1e-20; % keep old clusters, but assign them neglitable weight
      if(det(gmmf(i).sigma(:,:,iMix))<1e-20), gmmf(i).sigma(:,:,iMix)=gmmf(i).sigma(:,:,iMix)+0.00001*eye(3); end;
    end
  end

  gmmf(i).pi = gmmf(i).pi / sum( gmmf(i).pi );

end

for i=1:halfTimeWindow, gmmf(i)=gmmf(1+halfTimeWindow); end;
for i=(nframes-halfTimeWindow+1):nframes, gmmf(i)=gmmf(nframes-halfTimeWindow); end;

% Now Doing for BG
halfTimeWindow=floor(opts.bgModelTimeWindow/2);
kbWindow=zeros(nrows,ncols,opts.bgModelTimeWindow,'uint8');

% i = corr frame for which we are fitting the new gmm
for i=(1+halfTimeWindow):(nframes-halfTimeWindow),
  iStart=i-halfTimeWindow;
  iEnd=i+halfTimeWindow;
  % iStart and iEnd are the range of frames to which we fit our gmm's
  % kfWindow is the labels given by the current fg gaussian to the pixels in this window
  for j=iStart:iEnd,
    % j indexes the frames of this window, to fill in kfWindow
    if(isequal(curGmmb(j),curGmmb(i))),
      kbWindow(:,:,1+j-iStart)=kb(:,:,j);
    else,
      %computeIndices=find(T(:,:,j)=128|T(:,:,j)=64|T(:,:,j)=192);
      computeIndices=find(T(:,:,j)==128|T(:,:,j)==64); % No need to consider T=64, as it is hard background
      [xNotReq,kbWindow(:,:,1+j-iStart)]=vg_assignGmm_computeE(C(:,:,:,j),computeIndices,curGmmb(i));
      clear xNotReq;
    end
  end

  % So now we have kfWindow containing all the frames in the window with foreground gausssian labels

  for iMix=1:opts.nmixtures
    CiMix=[];
    for j=1:opts.bgModelTimeWindow
      iFrame=i-halfTimeWindow+(j-1);
      idxi = find( kbWindow(:,:,j) == iMix & A(:,:,iFrame)==0 & (T(:,:,iFrame)==128|T(:,:,iFrame)==64)  );
      CiMix = [CiMix; nema_impixel( C(:,:,:,iFrame), idxi )];
    end
    if size(CiMix,1)~=0
      gmmb(i).mu( :, iMix ) = mean( CiMix, 1 )';
      if(size(CiMix,1)<=2), 
        gmmb(i).sigma(:,:,iMix)=0.00001*eye(3);
      else
        gmmb(i).sigma( :, :, iMix ) = cov( CiMix );% + eye( 3 ) * 25/255/255; % add noise
      end;
      % If the cov-matrix is badly scaled, add some noise. It can be badly scaled for 2 reasons:
      % <= 2 members got assigned to this gauusian
      % the covariance is actually very low, because the pixels correspond to a uniform region
      if(det(gmmb(i).sigma(:,:,iMix))<1e-20), gmmb(i).sigma(:,:,iMix)=gmmb(i).sigma(:,:,iMix)+0.00001*eye(3); end;
      gmmb(i).pi( iMix ) = size(CiMix,1);
    else
      gmmb(i)=curGmmb(i);
      gmmb(i).pi( iMix ) = 1e-20; % keep old clusters, but assign them neglitable weight
      if(det(gmmb(i).sigma(:,:,iMix))<1e-20), gmmb(i).sigma(:,:,iMix)=gmmb(i).sigma(:,:,iMix)+0.00001*eye(3); end;
    end
  end
  gmmb(i).pi = gmmb(i).pi / sum( gmmb(i).pi );

end

for i=1:halfTimeWindow, gmmb(i)=gmmb(1+halfTimeWindow); end;
for i=(nframes-halfTimeWindow+1):nframes, gmmb(i)=gmmb(nframes-halfTimeWindow); end;

function beta = namg_beta3D( C )
%% Compute the precision of the pairwise Gaussian as a function of
%  the expected image gradient.
%
% Input parameters:
%   C - Input image.
%
% Output parameters:
%   beta - Gaussian precision for pairwise functions.

for i=1:size(C,4)
  beta(i) = 0.5 / nema_expected_image_gradient( C(:,:,:,i), 8 );
end

function betaTime = namg_betaTime3D( C )
%% Compute the precision of the pairwise Gaussian as a function of
%  the expected image gradient across time
%
% Input parameters:
%   C - Input image volume
%
% Output parameters:
%   betaTime - Gaussian precision for pairwise functions.

[nrows,ncols,ndims,nframes]=size(C);
betaTime=[];
for i=1:nframes-1
  C1=C(:,:,:,i); C2=C(:,:,:,i+1);
  sumOfDiff=nema_sum( (C1-C2).^2);
  avgTimeGrad=sumOfDiff/(nrows*ncols);
  betaTime(i) = 0.5 / avgTimeGrad;
end

function [A, energy, flow] = vg_dgraphcut( dgcHandle, STE, nrows,ncols,nframes)
%% Compute the minimum cut 
%
% Input Arguments
% dgcHandle -> handle to the graph to be passed to mexDGC
% STE -> the unary energies
%
% Output parameters:
%   A - The binary segmentation volume
%   energy - The energy of the segmentation [Currently disabled, just returns 0 here]
%   flow - The flow of the segmentation.

  %global scale,initInfo;
  global scale;
  [cut,flow] = mexDGC('minimizeDynamic',dgcHandle,uint32([1:size(STE,2)]),int32(round(STE*scale)),uint32([]),int32([]));
  A = reshape( ~double( cut ), nrows, ncols, nframes );
  %energy = namg_get_E( E_n, E_w * scale, STE * scale, ~A );
  energy=0;
