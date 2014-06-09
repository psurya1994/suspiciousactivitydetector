function [  ] = norm_angles2( angles )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

load dataset;

angles_n = angles-angles(1);
    
for i=1:5
    if(angles_n(i)<-180)
        angles_n(i) = angles_n(i) + 360;
    end
    if (angles_n(i)>=180)
        angles_n(i) = angles_n(i) - 360;
    end
end

for i=1:5
    if(angles(i)>=0)
        angles_mirror(i) = 180-angles(i);
    else
        angles_mirror(i) = -(180 + angles(i));
    end
end
angles_mirror_n = angles_mirror - angles_mirror(1);
angles_mirror_swap_n = [angles_mirror_n(1) angles_mirror_n(3) angles_mirror_n(2)...
                        angles_mirror_n(5) angles_mirror_n(4)];


ground_truth = [ground_truth;angles_n;angles_mirror_swap_n];

save ('dataset.mat', 'ground_truth');

end

