function [ suspicious ] = classify( angles, joints )

    suspicious = 0;
    if(angles(1)<80 || angles(1)>100)
        suspicious = 1;
    else
        if(~(angles(2)>90 && angles(2)<150)) % lifted hand, not in safe zone
            if((joints(1,4)<joints(1,2))&&(joints(2,4)<joints(2,2)))
                suspicious = 1; % break here
            else
                suspicious = 0;
            end
        else
            suspicious = 0;
        end
        if(~(angles(3)>30 && angles(3)<90)) % right hand, not in safe zone
            if((joints(1,3)<joints(1,5))&&(joints(2,3)>joints(2,5)))
                suspicious = 1; % break here
            else
                suspicious = (suspicious)||0;
            end
        else
            suspicious = (suspicious)||0;
        end
    end

end

