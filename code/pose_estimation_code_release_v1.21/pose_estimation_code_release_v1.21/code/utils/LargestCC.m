function CC = LargestCC(I)

% largest connected component of binary image I
%

[CCs n] = bwlabel(I);

best_cc = 0;
best_area = 0;
for ccix = 1:n
  area = length(find(CCs(:)==ccix));
  if area > best_area
    best_area = area;
    best_cc = ccix;
  end
end

CC = (CCs == best_cc);
