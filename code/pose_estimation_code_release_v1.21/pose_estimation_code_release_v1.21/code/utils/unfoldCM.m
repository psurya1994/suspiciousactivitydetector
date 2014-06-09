function CM = unfoldCM(CM,params)
% uncompreeses a given color model
% prevents color models usage with parameters different from those that were used to compute CM
% params.quantiz & params.colorSpace required
  if ~isfield(CM,'uniqueidx')
    return %already unfolded
  end
  if sum(params.quantiz == CM.quantiz) == 3 && strcmp(CM.colorSpace,params.colorSpace) % matching quantization and colorSpace with current params
    [cm Pfg Pbg] = deal(zeros(prod(CM.quantiz),size(CM.cm,2)));
    cm(CM.uniqueidx,:) = CM.cm;
    Pfg(CM.uniqueidx,:) = CM.Pfg;
    Pbg(CM.uniqueidx,:) = CM.Pbg;
    CM.cm = cm;
    CM.Pfg = Pfg;
    CM.Pbg = Pbg;
    CM = rmfield(CM,'uniqueidx');
  else
    error('cm have been computed with different set of parameters')
  end

end
  