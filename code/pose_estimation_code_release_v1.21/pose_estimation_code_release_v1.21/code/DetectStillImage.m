function angles = DetectStillImage(fullimgpath,pffubfmodel_path,facemodel_path,det_pars,verbose)
% runs upper body detector and optionally a opencv face detector on an image,
% face detection are regressed to the upper body detector coordinate frame using prelearned parameters
% Input:
% fullimgpath - relative/absolute path to an image
% ubfmodel_path - relative/absolute path to the pretrained upper body part-based model 
% facemodel_path - (optional) relative/absolute path to the pretrained opencv face model (xml file)
%                  if [] then skip face detection
% det_pars.iou_thresh - Intersection over Union threshold used during non-maximal suppression
%         .opencv_face_regparams - [p1 p2 p3 p4] - opencv_face detection to upperbody detection regression parameters
% verbose - 0 - no output
%         - 1 - print on screen
%         - 2 - show images
% Output:
% bbox(i,:)= [x1 y1 x2 y2] set of detections
%

  if nargin < 5
    verbose = 0;
  end
 
  if ~isempty(pffubfmodel_path)
    pffmodel = load(pffubfmodel_path); pffmodel = pffmodel.model;

    img = imread(fullimgpath);

    % find fine fixed aspect-ratio boxes;
    % upper body detection
    if isfield(det_pars,'ubfpff_scale') 
      tempimg = imresize(img,det_pars.ubfpff_scale);
    else
      tempimg = img;
    end

    boxes = detect(tempimg, pffmodel,det_pars.ubfpff_thresh);
    bbox = getboxes(pffmodel, boxes);

    if isfield(det_pars,'ubfpff_scale')
      bbox(:,1:end-1) = bbox(:,1:end-1)/det_pars.ubfpff_scale;
    end

    if ~isempty(bbox)
      det_hwratio = 0.9; %fixed aspect ratio
      bbox_center = (bbox(:,1:2)+bbox(:,3:4))/2;
      bbox_height = bbox(:,4)-bbox(:,2);
      dist_from_center = [1/det_hwratio*bbox_height/2 bbox_height/2];
      bbox = [bbox_center-dist_from_center bbox_center+dist_from_center bbox(:,end)];
      %bbox = boxes(:,[1:4 end]); % use the root boxes as the final prediction - aspect ratio is fixed v1.01 - less precise
      bbox = me_iou_nms(false, bbox, det_pars.iou_thresh);
    end
    if verbose 
      disp(['Image: ' fullimgpath]);
      disp(['   ' num2str(size(bbox,1)) ' upper bodies detected']);
    end
  else
    bbox = [];
  end
  
  % face_detection
  if ~isempty(facemodel_path)
    if exist('me_HaarDetectOpenCV') == 3 % if there exist mexed function
      fbox = me_HaarDetectOpenCV(facemodel_path,rgb2gray(gray2rgb(img)),30,30,1.1,det_pars.opencvface_minneighbours);% face detections in [x y w h] coor
      if ~isempty(fbox) % at least one face detected
        % regress faces to ubfdetections
        regp = det_pars.opencv_face_regparams;
        newx = fbox(:,1) - round(fbox(:,3) * regp(1));
        newy = fbox(:,2) - round(fbox(:,4) * regp(2));
        neww = fbox(:,3) * regp(3);
        newh = neww * regp(4);
        % convert to [x1 y1 x2 y2] and append a score -> use faces only as a complementary detector -> assign 10% lower score then min upper body score
        regubffromface = [newx newy newx+neww-1 newy+newh-1 repmat(det_pars.ubfpff_thresh-abs(det_pars.ubfpff_thresh/10),size(fbox,1),1)];
        % concatenate with ubdet & run nms
        bbox = [bbox; regubffromface];
        bbox = me_iou_nms(false, bbox, det_pars.iou_thresh);
      end
      if verbose 
        disp(['   ' num2str(size(fbox,1)) ' faces detected']);
      end
    else
      disp(['WARNING: opencv haardetection interface not found. Face model passed to the routine, but the face detetor is skiped.']);
    end
  end
  
  if verbose 
    disp(['   ' num2str(size(bbox,1)) ' detections kept after non-maximal suppression']);
  end
  
  if verbose > 1
    showboxes(img,bbox);
    axis equal;
    axis off;
  end
  
   %% Additional Code, manually added.
    
    disp('Number of detections = ');
    disp(size(bbox,1));
    startup;
    
    back = round(bbox(1:4)');
    back2 = [back(1) back(2) back(3)-back(1) back(4)-back(2)];
    
    [T sticks_imgcoor angles joints] = PoseEstimStillImage(pwd, 'images', '%06d.jpg', 1, 'ubf', back2', fghigh_params, parse_params_Buffy3and4andPascal, [], pm2segms_params, true);

  
end
