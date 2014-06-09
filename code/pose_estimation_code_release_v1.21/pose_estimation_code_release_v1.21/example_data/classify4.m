function [ similarity ] = classify4( angles, ground_truth )
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
        disp('angles_n = '); disp(angles_n);
        disp('ground_truth(i) = '); disp(ground_truth(i,:));
        deviation = 0;
        for j = 1:5
            dev = abs(angles_n(j)-ground_truth(i,j));
            if(dev>180)
                dev = 360-dev;
            end
            deviation = deviation + dev;
        end
        similarity = [similarity deviation];
    end
    
    disp(similarity);

end

