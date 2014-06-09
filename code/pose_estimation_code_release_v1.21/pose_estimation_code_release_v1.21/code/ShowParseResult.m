function fig = ShowParseResult(pm, show_whole, show_parts)

% draw posterior map results
% from Ramanan's human body parser
%

if nargin < 3
  show_parts = true;
end

if nargin < 2
  show_whole = true;
end

% if pm has been compressed, uncompress it
pm = UncompressPM(pm);

% show whole body
if show_whole
figure; subplot(1,2,1);
%pm.a = pm.a(:,:,[1 3 2]);                   % for paper figs (swap green/blue)
imagesc(uint8(pm.a*2e3));                   % keep org color range from Ramanan
axis equal; axis tight; axis off;
title('Whole body\newline pixel confidences color coded by part type');
%
subplot(1,2,2);
ppc_map = max(pm.b, [], 3);                 % per-pixel confidence map
imagesc(ppc_map, [0 1]);                    % color range reflect actual confidences
axis equal; axis tight; axis off;
title('Whole body\newline per-pixel confidences (max over parts)');
end

% show parts
% ids: 1 = torso; 2-5 = forearms and forelegs; 6-9 = end arms and end legs; 10 = head
if show_parts
Nparts = size(pm.b,3);
if Nparts == 10
  part_names = {'torso', 'upper left arm', 'upper right arm', 'upper left leg', 'upper right leg', ...
                'lower left arm', 'lower right arm', 'lower left leg', 'lower right leg', 'head'};
elseif Nparts == 6
 part_names = {'torso', 'upper left arm', 'upper right arm', 'lower left arm', 'lower right arm', 'head'};
else
  error('unknown model');
end
fig = figure;
nsf = ceil(sqrt(Nparts));
for spn = 1:Nparts            
  subplot(nsf,nsf,spn);
  imagesc(pm.b(:,:,spn), [0 1]);            % color range reflects actual pixel probs
  axis equal; axis tight; axis off;
  title(['part ' num2str(spn) '\newline(' part_names{spn} ')']);
end
set(fig, 'name', 'Body part segmentations'); 
drawnow;
end
