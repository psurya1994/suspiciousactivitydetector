function [ suspicious ] = classify( angles, joints )
% Input variables: angles- it is a vector containing angles of elevations
% of backbone, left upper limb, right upper limb, left lower limb and right
% lower limb in order (1x5 matrix).
% joints- Values of x and y coordinates of all joints ( 2x5 matrix).

% Output: returns '1' if the posture is suspicious.

   suspicious = 0;
    if(angles(1)<80 || angles(1)>100)  %if the person is bent more than 10 degrees from the upright position, it is suspicious.
        suspicious = 1;
    else
        if(~(angles(2)>90 && angles(2)<150)) % left arm not in safe zone;If the left arm is lifted more than 60 degrees from the backbone, it is suspicious.
            if((joints(1,4)<joints(1,2))&&(joints(2,4)<joints(2,2)))   % check orientation of fore arm
                suspicious = 1; % break here
            else
                suspicious = 0;
            end
        else
            suspicious = 0;
        end
        if(~(angles(3)>30 && angles(3)<90)) % right hand, not in safe zone; If the right arm is lifted more than 60 degrees from the backbone, it is suspicious.
            if((joints(1,3)<joints(1,5))&&(joints(2,3)>joints(2,5)))
                suspicious = 1; % break here
            else
                suspicious = (suspicious)||0;   %perform OR operation in order to prevent overwriting of the value.
            end
        else
            suspicious = (suspicious)||0;
        end
    end

end