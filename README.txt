
Steps to run the code.

1. Open MATLAB. Change the directory to the suspicious_activity_detector_v1 folder.

2. Add all the folders and subfolders in code directory to the MATLAB path.

3. Run the following commands one after another.
load ('detenv.mat');
startup;​​
cd code/voc-release3.1/voc-release3.1/
compile; % This doesn't work on windows. Suitable modification must be done.
cd ../../..

4. Change the threshold values for detection
det_pars.ubfpff_scale = 3;
det_pars.ubfpff_thresh = -0.75;
det_pars.iou_thresh = 0.9;

5. Code is setup and ready to run on a test image. Run the following commands to test.

image = imread('test_images/img_main.jpg');
[ubfdetections] = DetectStillImage2(image, 'pff_model_upperbody_final.mat', 'haarcascade_frontalface_alt2.xml', det_pars, 2);

6. If everything works fine, the result should be as shown in the attached image.

NOTE: 
The code currently works only on linux systems. Suitable modifications can be made for windows.

Debugging:
If you get the error that MATLAB is out of memory, re-run the code by reducing the value of det_pars.ubfpff_scale.
