function [detections Files] = DetectDir(root_dir,img_dir,pffubfmodel_path,facemodel_path,det_pars,draw,verbose)
% executes DetectStillImage on a directory with images
% Input:
%     img_dir - relative/absolute path to an directory with images
%     pffubfmodel_path - relative/absolute path to the pretrained upper body part-based model 
%     facemodel_path - (optional) relative/absolute path to the pretrained opencv face model (xml file)
%                      if [] then skip face detection
%     det_pars.iou_thresh - Intersection over Union threshold used during non-maximal suppression
%             .opencv_face_regparams - [p1 p2 p3 p4] - opencv_face detection to upperbody detection regression parameters
%     verbose - 0 - no output
%                - 1 - print on screen    
%                - 2 - show images
%  Output:
%     detections - detection results for a directory, detections{i} - detections found in Files(i).name
%     Files - struct array as return by dir matlab function for all processed images.
%  routine creates also a subdirectory 'dets' where results are saved
    
    Files = dir(fullfile(root_dir,img_dir));
    invalid = false(length(Files),1);
    RegularExpression = '(\w+\.(jpg)|(jpeg)|(gif)|(bmp)|(png)|(ppm))$';
    for i=1:numel(Files)
      invalid(i) = isempty(regexpi(Files(i).name, RegularExpression));
    end
    Files(invalid) = [];
    
    detections = cell(1,length(Files));
    
    dets_dir = 'calvinubf_dets';
    
    if ~exist(fullfile(root_dir,dets_dir),'dir')
      mkdir(fullfile(root_dir,dets_dir));
    end

    idxes = 1:numel(Files);
    
    for idx=idxes
      if verbose > 1
        clf;
      end
      detections{idx} = DetectStillImage(fullfile(root_dir,img_dir,Files(idx).name),pffubfmodel_path,facemodel_path,det_pars,verbose);   
      if verbose > 1
        ginput(1);
      end
      if draw
        img = imread(fullfile(root_dir,img_dir,Files(idx).name));
        bb = detections{idx};
        if ~isempty(bb)
           bb(:,3:4) = bb(:,3:4) - bb(:,1:2) +1; %convert to [x y w h] for drawing
           for d=1:size(bb,1)
             img = PaintBB(img,round(bb(d,1:4)),[1 0 0],[1 2 3]);
           end
        end
        [trash fname fext] = fileparts(Files(idx).name);
        detfile = fullfile(root_dir,dets_dir,[fname '_d' fext]);
        imwrite(img,detfile);
      end
    end
  
end