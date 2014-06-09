function [ angles joints ] = extract_stick( curSticks )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
%% FInding exact angles 
  jointsX = [];
  jointsY = [];
  % body
  body = (curSticks(:,:,1)==255).*(curSticks(:,:,2)==0).*(curSticks(:,:,3)==0);
  k =  my_hough(double(body),2);
  theta = k(1).theta+90;
  
  if(k(1).point1(2)>k(1).point2(2))
      jointsX = [jointsX k(1).point2(1)];
      jointsY = [jointsY k(1).point2(2)];
  else
      jointsX = [jointsX k(1).point1(1)];
      jointsY = [jointsY k(1).point1(2)];
  end
  

  %upper limbs
  ulimb = (curSticks(:,:,1)==0).*(curSticks(:,:,2)==255).*(curSticks(:,:,3)==0);
  k1 = my_hough(ulimb,8);
  pointX=[];
  for i=1:8
    pointX = [pointX k1(i).point1(1)];
  end
  
  avgX=mean(pointX);

    for i=1:8
        if(pointX(i)<avgX)
            lultheta = k1(i).theta+90;
            if(distance(jointsX(1),jointsY(1),k1(i).point1(1),k1(i).point1(2))>distance(jointsX(1),jointsY(1),k1(i).point2(1),k1(i).point2(2)))
                jointsX = [jointsX k1(i).point1(1)];
                jointsY = [jointsY k1(i).point1(2)];
            else
                jointsX = [jointsX k1(i).point2(1)];
                jointsY = [jointsY k1(i).point2(2)];
            end
            break;
        end;
    end;
  
  
    for i=1:8
        if(pointX(i)>avgX)
            rultheta = k1(i).theta+90;
            if(distance(jointsX(1),jointsY(1),k1(i).point1(1),k1(i).point1(2))>distance(jointsX(1),jointsY(1),k1(i).point2(1),k1(i).point2(2)))
                jointsX = [jointsX k1(i).point1(1)];
                jointsY = [jointsY k1(i).point1(2)];
            else
                jointsX = [jointsX k1(i).point2(1)];
                jointsY = [jointsY k1(i).point2(2)];
            end
            break;
        end;
    end;

    
    
    % lower limbs
    
    llimb = (curSticks(:,:,1)==255).*(curSticks(:,:,2)==255).*(curSticks(:,:,3)==0);
    
    k2 = my_hough(llimb,8);
    pointY=[];
    for i=1:8
        pointY = [pointY k2(i).point1(1)];
    end
  
    avgY=mean(pointY);

    
    
    for i=1:8
        if(pointY(i)<avgY)
            llltheta = k2(i).theta+90;
            if(distance(jointsX(2),jointsY(2),k2(i).point1(1),k2(i).point1(2))>distance(jointsX(2),jointsY(2),k2(i).point2(1),k2(i).point2(2)))
                jointsX = [jointsX k2(i).point1(1)];
                jointsY = [jointsY k2(i).point1(2)];
            else
                jointsX = [jointsX k2(i).point2(1)];
                jointsY = [jointsY k2(i).point2(2)];
            end
            break;
        end;
    end;
    
    
    for i=1:8
        if(pointY(i)>avgY)
            rlltheta = k2(i).theta+90;
            if(distance(jointsX(3),jointsY(3),k2(i).point1(1),k2(i).point1(2))>distance(jointsX(3),jointsY(3),k2(i).point2(1),k2(i).point2(2)))
                jointsX = [jointsX k2(i).point1(1)];
                jointsY = [jointsY k2(i).point1(2)];
            else
                jointsX = [jointsX k2(i).point2(1)];
                jointsY = [jointsY k2(i).point2(2)];
            end
            break;
        end;
    end;

  
    angles = [theta lultheta rultheta llltheta rlltheta];
    angles_n = angles - theta;
    joints = [jointsX;jointsY];

end

