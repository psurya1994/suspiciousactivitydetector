function [MAP, segsAll, respIm] = sampleFromRespIm(respIm, model, nsamples, prob_exp_smoothing)
% sampling from joint distribution 
% for tree-structured models in model
% realized by sampling at a root marginal and then sampling at conditional distributions of the children.

% size(respIm) must be 4 dimensional
%
% model.orient must all be in [0,1] and sum-to-1
%
% Adaptation of Ramanan's code by M.Eichner
%
% !! Important: respIm must play properly with model.mu !
% In Ramanan's models, model.mu assuems respIm has segments
% anchored on middle-top (and not center)
%
% nsamples - number of samples to be generated 
% prob-exp-smoothing -  if < 1 smoothing the probability at every stage before sampling 
%                       if > 1 sharpening 
%
% returned data-structure is a cell array of size nsamples+1 with parameters describing sticks
% the first cell array entry is the 'BEST' pose - (at each level instead of sampling from the distribution we select argmax)

unary = respIm;

if nargin < 3 
  nsamples = 0;
end

if nargin < 4 || isempty(prob_exp_smoothing)
  prob_exp_smoothing = 1;
end

% returns nsamples samples from the joint where fist sample is 'MAP' 

% shortcuts
[imy imx trash trash] = size(respIm);              % keep trash, as otherwise imy imx corrupted
numClass = max(model.equiv_class);
numTypes = length(model.parents);                  % number of body parts
imo = size(model.orient,2);                        % number of possible limb orientations
assert(numClass == numTypes);
assert(length(size(respIm))==4);
assert(all(model.orient(:)>=0));

% message passing optimization
% updates all limbs respIm
% by optimizing energy fct of
% unary terms (image evidence above) and
% pairwise terms (rel loc/orient betw connected limbs, with pre-learnt priors)
% output is updated respIm(x,y,orient,limb_type)
mean_y = round(model.mu(:,2));
mean_x = round(model.mu(:,1));
for p = model.post(1:end-1),
  parType = model.parents{p};
  
  currResp = respIm(:,:,:,p);
  
  currResp = reshape(currResp,imy*imx,imo);

  % Reflect borders (number of orientations becomes 2*imo-1)
  currResp = cat(2,currResp(:,(imo/2+1):end),currResp,currResp(:,1:(imo/2-1)));

  % Convolve (apply relative orientation part of kinem prior)
  currResp = conv2(currResp,fliplr(model.orient(p,:)),'valid');
  currResp = reshape(currResp,[imy imx imo]);  
    
  % Accumulate support for the parent)
  currResp = local_sum_zero(currResp,model.box(p)*2);
  
  % Shift resp to the position of the parent anchor
  currResp = partshiftZ0(currResp,mean_y(p),mean_x(p));
  
  % Pass message to the parent
  respIm(:,:,:,parType) = respIm(:,:,:,parType).*currResp;
end

%at the root
p = model.root;
currResp = respIm(:,:,:,p);

% absolute orientation prior is applied already to unaries don't need to apply it second time 
%currResp = currResp.*repmat(permute(model.orient(p,:),[1 3 2]),[imy imx 1]);

% renormalize
NN = sum(currResp(:));
RootMarginal = currResp/NN;


% segsAll = struct('u',{},'v',{},'len',{},'wid',{},'x',{},'y',{},'o',{},...
%                              'type',{},'prob',{},'sticks',{});
segsAll = struct('x',{},'y',{},'o',{},'prob',{},'sticks',{});
%segsAll = struct('o',{},'prob',{},'sticks',{});


%sampling procedure - start by sampling from the marginal of the root and then sample from conditionals


for big_iter = 1:nsamples + 1,
  
    p = model.root;
  
    [Y,X,I,P] = deal(zeros(numTypes,1));
  
    if big_iter == 1,
        BEST = 1;
    else
        BEST = 0;
    end

    if BEST,
        [trash,ii] = max(RootMarginal(:));
    else
        %Sample
        ii = find(sampleWithR(RootMarginal(:).^prob_exp_smoothing,2),1);
    end
    [Y(p) X(p) I(p)] = ind2sub(size(RootMarginal),ii);
    P(p) = RootMarginal(ii);

        % Pass messages down
    % Define downstream orientations as walking counter-clockwise starting 
    % from 13, pointing down
    % ie [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24]
    % -> [1 24 23 22 21 20 19 18 17 16 15 14 13 12 11 10 9 8 7 6 5 4 3 2]
    % this routine is identical to that in the first parse
    for p = model.pre(2:end),
        parType = model.parents{p};
        % as orient prior may not be symmetric then turn it around for messages from parent to child
        currOrient = model.orient(p,[1 end:-1:2]);
           
        currCoor = coordshiftZ0([Y(parType),X(parType),I(parType)],size(respIm),-mean_y(p),-mean_x(p));
        subcube = zeros(1, 1,imo);
        subcube(1,1,currCoor(3)) = 1;
        subcube = reshape(subcube,1,imo);
        subcube = cat(2,subcube(:,(imo/2+1):end),subcube,subcube(:,1:(imo/2-1)));
        
        subcube = conv2(subcube,fliplr(currOrient),'valid');
        subcube = reshape(subcube,[1 1 imo]);
        
        sc_y1 = max(currCoor(1)-model.box(p), 1);
        sc_y2 = min(currCoor(1)+model.box(p)-1, imy);
        sc_x1 = max(currCoor(2)-model.box(p), 1);
        sc_x2 = min(currCoor(2)+model.box(p)-1, imx);
        
        subcube = repmat(subcube,[sc_y2-sc_y1+1 sc_x2-sc_x1+1]);
        % each node (respIm) is already multiplied by all the incoming messages from its children
        % so we can just multiply it by conditional from its parent
        subcube = respIm(sc_y1:sc_y2, sc_x1:sc_x2, :, p) .* subcube;

        % by the definition the states outside the cube have 0 probability
        subcube = subcube / sum(subcube(:)); %renormalize the conditional probability to sum to one
        if BEST,
            [trash,ii] = max(subcube(:));
        else
           %Sample
           ii = find(sampleWithR(subcube(:).^prob_exp_smoothing,2),1);
        end
        [offset(1) offset(2) offset(3)] = ind2sub(size(subcube),ii);
        Y(p) = sc_y1+offset(1)-1;
        X(p) = sc_x1+offset(2)-1;
        I(p) = offset(3);
        P(p) = subcube(ii);
   
    end
    
    % calc prob of the configuration
    prob = prod(P);
    
    %Shift the states to align with center of segment
    for p = 1:numTypes,
      currCoor = coordshiftZ0([Y(p) X(p) I(p)],size(respIm),-model.len(p), 0);
      Y(p) = currCoor(1);
      X(p) = currCoor(2);
      I(p) = currCoor(3);
    end

    %Convert to segments
    stepo = 360/imo;
    [uu,vv] = pol2cart(([stepo:stepo:360] + 90 - stepo + 180) * pi/180,1);
    vv = -vv;

    U = uu(I)';
    V = vv(I)';
    Len = model.len;
    p1 = [X Y] - [U V].*repmat(Len,1,2); % body anchor point
    p2 = [X Y] + [U V].*repmat(Len,1,2); 

    if BEST
      MAP.x = single(X); % center x
      MAP.y = single(Y); % center y
      MAP.o = single(I); % orientation
      MAP.prob = prob; % probability of configuration 
      MAP.sticks = single(([p1'; p2'])); % end points of segments representing body parts
    else
      %seg.type = [1:numTypes]'; 
      seg.x = single(X); % center x
      seg.y = single(Y); % center y
      seg.o = single(I); % orientation
      seg.prob = prob; % probability of configuration 
      seg.sticks = single(([p1'; p2'])); % end points of segments representing body parts
      segsAll(big_iter-1) = seg;
    end
    
end


end

