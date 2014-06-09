function CM = foldCM(CM,params)
% compresses the color model, 
% make sure that the CM is complete
% params.quantiz and params colorSpace required



  uniqueidx =  find(sum(abs(CM.Pbg)+abs(CM.Pfg),2) ~= 0); % find all bin indexes where fg or bg models are nonzero
  CM.cm=CM.cm(uniqueidx,:);
  CM.Pfg=CM.Pfg(uniqueidx,:);
  CM.Pbg=CM.Pbg(uniqueidx,:);
  CM.uniqueidx = uniqueidx;
  CM.quantiz = params.quantiz;
  CM.colorSpace = params.colorSpace;
  end
  