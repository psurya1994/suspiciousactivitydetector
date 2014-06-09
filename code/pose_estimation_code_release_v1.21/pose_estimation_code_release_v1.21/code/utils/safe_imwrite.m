function safe_imwrite(img, imgname)
max_iter = 4;
timedelay = 20;
ok = false;
attempt = 1;
while ~ok && attempt <= max_iter
  try 
    imwrite(img,imgname);
    ok = true;
  catch
    attempt = attempt + 1;
    pause(timedelay);
  end
end

if attempt > max_iter
  disp(['WARNING: image ' imgname ' could not be saved']);
end

