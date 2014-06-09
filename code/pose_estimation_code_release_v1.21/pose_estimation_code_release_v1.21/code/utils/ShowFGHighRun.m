function fid = ShowFGHighRun(STD, STM, p)

% displays a the video sequence STD(:,:,:,t) and
% corresponding labelling STM(:,:,t)
%
% pause for p seconds between frames (p can be fractional)
% return figure handler fid with videos
%
% if p<0 -> wait for user return
%

fid = figure;
for cur_dix = 1:size(STD,4)
    subplot(1,2,1); vimagesc(STD(:,:,:,cur_dix));
    subplot(1,2,2); vimagesc(STM(:,:,cur_dix));
    set(fid,'name',['run to be segmented - ' num2str(cur_dix) 'th frame in the run']);
    drawnow;
    if p >= 0
      pause(p);
    else
      keyboard;
    end
end
