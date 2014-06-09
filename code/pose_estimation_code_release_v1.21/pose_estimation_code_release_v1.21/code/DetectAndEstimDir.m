function DetectAndEstimDir(img_dir,pffubfmodel_path,facemodel_path,det_pars,classname,fghigh_pars,parse_pars,addinf_pars,segm_pars,verbose)
    Files = dir(img_dir);
    invalid = false(length(Files),1);
    RegularExpression = '(\w+\.(jpg)|(jpeg)|(gif)|(bmp)|(png)|(ppm))$';
    for i=1:numel(Files)
      invalid(i) = isempty(regexpi(Files(i).name, RegularExpression));
    end
    Files(invalid) = [];
    
    %parse_pars.use_fg_high = false % uncomment this line if you want to skip the foreground highlighting stage
 
    for idx=1:numel(Files);
      
      stick_coor = cell(0);
      T = struct('D',{}, 'FGH',{}, 'PM',{},'CM',{});

      
      [trash,imgname,imgext] = fileparts(Files(idx).name);
      outname_mat = fullfile(img_dir,[imgname imgext '_' classname '_pms.mat']);
      if exist(outname_mat,'file')
        continue;
      end
      
      detections = DetectStillImage(fullfile(img_dir,Files(idx).name),pffubfmodel_path,facemodel_path,det_pars,verbose);
      
      if ~isempty(detections)
        %convert coordinates from [x y width height] to [x1 y1 x2 y2]
        detections(:,3:4) = detections(:,3:4) - detections(:,1:2) +1;

        img_name_format = [imgname '_%d' imgext];
        temp_dir = 'temp';
        if ~exist(fullfile(img_dir,temp_dir),'dir')
          mkdir(fullfile(img_dir,temp_dir));
        end
        dets_dir = 'dets';
        if ~exist(fullfile(img_dir,dets_dir),'dir')
          mkdir(fullfile(img_dir,dets_dir));
        end

        img = imread(fullfile(img_dir,Files(idx).name));

        for dix=1:size(detections,1) % run pose estimation for every detection
          imwrite(img,fullfile(img_dir,temp_dir,sprintf(img_name_format,dix)));
          [T(dix) stick_coor{dix}] = PoseEstimStillImage(img_dir,temp_dir,img_name_format,dix, classname, round(detections(dix,1:4)'), fghigh_pars, parse_pars, addinf_pars, segm_pars, verbose);
          delete(fullfile(img_dir,temp_dir,sprintf(img_name_format,dix)));
        end
        
        if ~isempty(detections)
          for d=1:size(detections,1)
            img = PaintBB(img,round(detections(d,1:4)),[1 0 0],[1 2 3]);
          end
        end
        detfile = fullfile(img_dir,dets_dir,[imgname '_d' imgext]);
        imwrite(img,detfile);
      end
      save(outname_mat,'T','stick_coor','detections');
    end

    
end