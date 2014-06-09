function intimg = integralImageXY(img)
  intimg = cumsum(cumsum(img,1),2);
  