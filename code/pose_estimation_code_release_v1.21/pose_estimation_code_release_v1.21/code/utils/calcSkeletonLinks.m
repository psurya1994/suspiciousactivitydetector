function [links sticks] = calcSkeletonLinks(sticks,orderedsticks,rootidx,top_rootidx,topchildren,bottomchildren)
% finds links between body parts 
% and try to sort sticks so that sticks(1:2,i) is the anchor point of the stick wrt to its parent
%    sticks(:,N) - matrix of N sticks, 
%               where sticks(:,i) = [x1 y1 x2 y2]' are the coordinates of endpoints of i-th stick
%    orderedsticks - true -> follow the tree as sticks(1:2,:) is always an anchor point
%                  - false -> figure out which is the anchor point as traversing through the tree 
%    rootidx - root index in sticks
%    top_rootidx - coor index {1,2} from which topchildren tree is traversed 
%                  the remaining coor is then assigned to bottom_rootidx 
%    topchildren bottomchildren are the trees traversed from top_rootidx and bottom_rootidx respectively

  Nsticks = size(sticks,2);
  if nargin < 3
    topchildren = cell(1,Nsticks);
    bottomchildren = cell(1,Nsticks);
    switch Nsticks
      case 6
        topchildren{1} = [2 3 6];
        topchildren{2} = 4;
        topchildren{3} = 5;
      case 10
        topchildren{1} = [2 3 10];
        topchildren{2} = 6;
        topchildren{3} = 7;
        bottomchildren{1} = [4 5];
        bottomchildren{4} = 8;
        bottomchildren{5} = 9;
      otherwise
        error('unsupported class type')
    end
    rootidx = 1;
    [trash top_rootidx] = min(sticks([2 4],rootidx)); 
  end

  links = cell(0);
  
  if ~orderedsticks % sticks are not ordered

    bottom_rootidx = mod(top_rootidx,2)+1;
    usedcoor = zeros(2,Nsticks);


    if me_isEmptyStick(sticks(:,rootidx)) % if root is occluded then no links between body parts can be established
      links = [];
      return
    end

    idx = 0;
    [links usedcoor idx] = findlink(idx,links,usedcoor,sticks,topchildren,rootidx,top_rootidx);
    [links usedcoor idx] = findlink(idx,links,usedcoor,sticks,bottomchildren,rootidx,bottom_rootidx);

    % flip coordinates so that the limb anchor 
    % wrt its parent is always the first coordinate
    for s=find(usedcoor(2,:)-usedcoor(1,:) < 0)
      sticks(:,s) = sticks([3 4 1 2],s);
    end

  else % sticks are ordered: always the top one is the anchor point
    
    occluded = me_isEmptyStick(sticks);
    
    % first create links for root's immediate topchildren
    for ch=1:numel(topchildren{rootidx});
      if occluded(ch)
        continue;
      end
      links{end+1} = [sticks(1:2,rootidx); sticks(1:2,topchildren{rootidx}(ch))];  
    end
    % create links for the root's immediate bottomchildren
    for ch=1:numel(bottomchildren{rootidx});
      if occluded(ch)
        continue;
      end
      links{end+1} = [sticks(3:4,rootidx); sticks(1:2,bottomchildren{rootidx}(ch))];  
    end
    % create all other links
    nonrootchildren = cellfun(@(x,y)[x y]...
      ,[topchildren([1:rootidx-1]) {[]} topchildren([rootidx+1:end])]...
      ,[bottomchildren([1:rootidx-1]) {[]}  bottomchildren([rootidx+1:end])]...
      ,'UniformOutput',false);
    for p=find(~cellfun(@isempty,nonrootchildren))
      if occluded(p)
        continue;
      end
      for ch=1:numel(nonrootchildren{p});
        if occluded(ch)
          continue;
        end
        links{end+1} = [sticks(3:4,p); sticks(1:2,nonrootchildren{p}(ch))];  
      end
    end
    % sticks already ordered no need to reorder them
  end
  
  links = [links{:}];

end


function [links usedcoor idx] = findlink(idx,links,usedcoor,sticks,tree,pix,cix)
  idx = idx+1;
   usedcoor(cix,pix) = idx;
   parent_start = iSwitch(mod(cix,2),1:2,3:4);
   for child=tree{pix}
     if ~me_isEmptyStick(sticks(:,child))
       cend = iSwitch(sum((sticks(parent_start,pix)-sticks(1:2,child)).^2) > sum((sticks(parent_start,pix)-sticks(3:4,child)).^2),2,1);
       idx = idx + 1;
       usedcoor(cend,child) = idx;
       child_end = iSwitch(mod(cend,2),1:2,3:4);
       links{end+1} = [sticks(parent_start,pix); sticks(child_end,child)];
       [links usedcoor idx] = findlink(idx,links,usedcoor,sticks,tree,child,mod(cend,2)+1);
     end
   end
end