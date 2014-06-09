function Marginal = quick_msgpass_tree_corr(respIm, model, shift, motion_shift)
% this is a corrected version of quick_msgpass_tree, 
% 'in place' message passing has been replaced with storing all the messages and at the end computing marginals 
% this produces correct marginals for all nodes (in quick_msgpass_tree correct was only marginal for the root node.)

% Rapid sum-product algorithm
% for tree-structured models in model
%
% size(respIm) must be length 4
%
% model.orient must all be in [0,1] and sum-to-1
%
% Adaptation of Ramanan's code by V. Ferrari
%
% shift: whether to shift the respIm to align with center of segment
% might be useful not to do it if multiple infer
%
% !! Important: respIm must play properly with model.mu !
% In Ramanan's models, model.mu assuems respIm has segments
% anchored on middle-top (and not center)
%
% !!! motion_shift should be only past in case of temporal inference across several frames for a single limb

% process arguments
if nargin < 4
  motion_shift = [];
end
if nargin < 3
  shift = true;
end

% shortcuts
[imy imx trash trash] = size(respIm);              % keep trash, as otherwise imy imx corrupted
numClass = max(model.equiv_class);
numTypes = length(model.parents);                  % number of body parts
imo = size(model.orient,2);                        % number of possible limb orientations
assert(numClass == numTypes);
assert(length(size(respIm))==4);
assert(all(model.orient(:)>=0));

incMsg = cell(numTypes); % keeps all messages incoming to nodes
Marginal = zeros(size(respIm));

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
  % Multiply in all the messages from kids to the current respIm
  % select current unary potential
  currResp = respIm(:,:,:,p);
  for i=1:numTypes
    if i ~= parType && ~isempty(incMsg{i,p})
      currResp = currResp .* incMsg{i,p};
    end
  end

   % if motion_shift is passed then shift the msg accordingly
  if ~isempty(motion_shift)
    % if motion shift is passed then the model is expected to contain a single body part in different time instances
    % and the structure of the graph is always a chain
    if parType < p 
      currResp = shiftRespYX(currResp,-motion_shift.Y(:,:,:,p),-motion_shift.X(:,:,:,p),true);
    else 
      currResp = shiftRespYX(currResp,motion_shift.Y(:,:,:,p),motion_shift.X(:,:,:,p),true);
    end
    %currResp = currResp/sum(currResp(:)); %when message is considered we should not renormalize it
  end
  
  currResp = reshape(currResp,imy*imx,imo);

  % Reflect borders (number of orientations becomes 2*imo-1)
  currResp = cat(2,currResp(:,(imo/2+1):end),currResp,currResp(:,1:(imo/2-1)));

  % Convolve (apply relative orientation part of kinem prior)
  currResp = conv2(currResp,fliplr(model.orient(p,:)),'valid');
  currResp = reshape(currResp,[imy imx imo]);  
    
  % Accumulate support for the parent)
  currResp = local_sum_zero(currResp,model.box(p)*2);
  
  % Shift resp to the position of the parent middle-top anchor
  currResp = partshiftZ0(currResp,mean_y(p),mean_x(p));
  
  % Pass message to the parent
  incMsg{p,parType} = currResp;
  
end

%at the root
p = model.root;
currResp = respIm(:,:,:,p);
%apply absolute orientation prior to the root
currResp = currResp.*repmat(permute(model.orient(p,:),[1 3 2]),[imy imx 1]);
% renormalize and store in the root likelihood
NN = sum(currResp(:));
respIm(:,:,:,p) = currResp/NN;

% Pass messages down
% Define downstream orientations as walking counter-clockwise starting 
% from 13, pointing down
% ie [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24]
% -> [1 24 23 22 21 20 19 18 17 16 15 14 13 12 11 10 9 8 7 6 5 4 3 2]
% this routine is identical to that in the first parse
for p = model.pre(2:end),
    % prepare the incoming message from the parent
    parType = model.parents{p};
    currResp = respIm(:,:,:,parType);
    % point-wise multiply the unary of the parent with all messeges incoming to the parent 
    % except the one from current branch
    for i=1:numTypes
      if i ~= p && ~isempty(incMsg{i,parType}) 
        currResp = currResp .* incMsg{i,parType};
      end
    end
    
    %Shift resp to the position in the childrens middle-top anchor
    currResp = partshiftZ0(currResp,-mean_y(p),-mean_x(p));
    
     % if motion_shift is passed then shift the msg accordingly
    if ~isempty(motion_shift)
      if parType < p
        currResp = shiftRespYX(currResp,motion_shift.Y(:,:,:,p),motion_shift.X(:,:,:,p),true);
      else 
        currResp = shiftRespYX(currResp,-motion_shift.Y(:,:,:,p),-motion_shift.X(:,:,:,p),true);
      end
    end
    
    % apply kinematic prior
    currResp = reshape(currResp,imy*imx,imo);

    % Now do normal message passing
    % Circularly convolve
    currResp = cat(2,currResp(:,(imo/2+1):end),currResp,currResp(:,1:(imo/2-1)));
    
    % as orient prior may not be symmetric then turn it around for messages from parent to child
    currOrient = model.orient(p,[1 end:-1:2]);
    currResp = conv2(currResp,fliplr(currOrient),'valid');
    currResp = reshape(currResp,[imy imx imo]);
    currResp = local_sum_zero(currResp,model.box(p)*2);
       
    
    % store the incoming message
    incMsg{parType,p} = currResp;
end

%calculate marginals
for p=1:numTypes
  % it means the unary of a node point-wise multiplied with all incoming messages
  currResp = respIm(:,:,:,p);
  for i=1:numTypes
    if ~isempty(incMsg{i,p})
      currResp = currResp .* incMsg{i,p};
    end
  end
  Marginal(:,:,:,p) = currResp/sum(currResp(:));
end


% Shift the respIm to align with center of segment
if shift
  Marginal = align_respIm(Marginal, model);
end
