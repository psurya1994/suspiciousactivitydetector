function DrawBB(BB, col, bb_thick,  score, id,  pars,  limits,  logo)
% DrawBB(BB, col, bb_thick,  score, id,  pars,  limits,  logo)
% Adds bounding-box BB to current figure
%
% BB in form:
% [minx maxx miny maxy]
%
% if scalar score given -> draw score in color score_col
% same for id; pars is a drawing parameters structure, necessary if wanting
% to draw more than plainly a BB
%
% if limits = [width height] given -> don't draw anything outside the image
% (good for exporting to files properly, when automatically cropping border)
%
% A logo can be drawn in the BB, pars should be defined
%     .logo_col
%     .logo_size
%
% MOD: mjmairn
%

%global colors;
colors = sampleColors;

if nargin < 3
  bb_thick = 1;
end

% resolves colors (which can be input as either true color 3x1 vectors,
% or indeces into colors(col_id,:)
if nargin < 2
  col = [0 1 0];            % green is the default color
elseif length(col)==1       % index in colors
  col_id = col;
  col_id = mod(col_id-1,size(colors,1))+1;
  col = colors(col_id,:);
end

if nargin >= 6              % settle color of score
if length(pars.score_col)==1   % index in colors
  col_id = pars.score_col;
  col_id = mod(col_id-1,size(colors,1))+1;
  pars.score_col = colors(col_id,:);
end
end

if nargin >= 6              % setttle color of id
if length(pars.id_col)==1   % index in colors
  col_id = pars.id_col;
  col_id = mod(col_id-1,size(colors,1))+1;
  pars.id_col = colors(col_id,:);
end
end

if nargin >= 6              % settle color of logo
if length(pars.logo_col)==1   % index in colors
  col_id = pars.logo_col;
  col_id = mod(col_id-1,size(colors,1))+1;
  pars.logo_col = colors(col_id,:);
end
end

if nargin < 7
  limits = false;
end


% Draw BB (= 4 line segments)
hold on;
if limits
  PlotBBClipped(BB, col, bb_thick, limits);
else
  plot([BB(1) BB(2) BB(2) BB(1) BB(1)], [BB(3) BB(3) BB(4) BB(4) BB(3)], 'Color', col, 'LineWidth', bb_thick);
end

if nargin < 4
  return;
end

% Draw score
if pars.draw_bb_score
  min_x = BB(1)+bb_thick/2+4;
  min_y = BB(3)+bb_thick/2+pars.score_size/2+3;
  % text(x,y,...) is drawn with the bottom-right of the first char on x,y
  if limits
    ok = (min_x > 0) & (min_y-pars.score_size > 0) & (min_x+pars.id_size*(floor(log10(id))+1) < limits(1)) & (min_y+pars.id_size < limits(2));
  else
    ok = true;
  end
  if ok
    text(min_x, min_y, num2str(score,'%1.3f'), 'color', pars.score_col, 'FontSize', pars.score_size, 'FontWeight', 'bold');
  end
end

% Draw id
if pars.draw_id
  min_x = BB(1)+bb_thick/2+4;
  min_y = BB(4)-bb_thick/2-pars.id_size/2;
  if limits
    ok = (min_x > 0) & (min_y-pars.score_size > 0) & (min_x+pars.id_size*(floor(log10(id))+1) < limits(1)) & (min_y+pars.id_size < limits(2));
  else
    ok = true;
  end
  if ok
    text(min_x, min_y, num2str(id), 'color', pars.id_col, 'FontSize', pars.id_size, 'FontWeight', 'bold');
  end
end

% Draw logo
if pars.draw_logo
  min_x = BB(2)-bb_thick-14;
  min_y = BB(3)+bb_thick/2+pars.score_size/2+3;
  if limits
    ok = (min_x > 0) & (min_y-pars.score_size > 0) & (min_x+pars.id_size*(floor(log10(id))+1) < limits(1)) & (min_y+pars.id_size < limits(2));
  else
    ok = true;
  end
  if ok
    text(min_x, min_y, logo, 'color', pars.logo_col, 'FontSize', pars.logo_size, 'FontWeight', 'bold');
  end
end
