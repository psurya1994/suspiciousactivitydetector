function cols = class_id2cols(id)

% colors for drawing limbs depending on the object class 
%

switch id
    case 1  % ubp
         cols = [];       % undefined yet
    case 2  % ubf
         %cols = [1 0 0; 0 1 0; 0 1 0; 0 0 1; 0 0 1; 1 0 1];  % on screen
         %cols = [1 1 0; 0 1 0; 0 1 0; 0.5 0.5 1; 0.5 0.5 1; 1 0 1];      % brighter -> for paper
         cols = [1 0 0; 0 1 0; 0 1 0; 1 1 0; 1 1 0; 1 0 1];  % paper
    case 3  % full
         cols = [1 0 0; 0 1 0; 0 1 0; 0 0 1; 0 0 1; 1 1 0; 1 1 0; 0 1 1; 0 1 1; 1 0 1];  % paper
    otherwise
        error([mfilename ': unknown class id ' num2str(id)]);
end
