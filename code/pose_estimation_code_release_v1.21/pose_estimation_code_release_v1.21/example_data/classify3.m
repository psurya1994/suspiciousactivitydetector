function [ similarity ] = classify3( angles, ground_truth )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
    angles_n = angles-angles(1);
    
    for i=1:5
        if(angles_n(i)<-180)
            angles_n(i) = angles_n(i) + 360;
        end
        if (angles_n(i)>=180)
            angles_n(i) = angles_n(i) - 360;
        end
    end
    similarity = [];
    for i=1:size(ground_truth,1)
        similarity = [similarity pdist2(angles_n,ground_truth(i,:))];
    end
    
    disp(similarity);

end

