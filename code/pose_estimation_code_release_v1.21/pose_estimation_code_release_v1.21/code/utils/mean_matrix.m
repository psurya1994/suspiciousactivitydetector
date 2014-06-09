function o = mean_matrix(M)
  
% matrix mean
% like mean,
% but returns mean over columns
% always, even if M is a vector

if size(M,1) == 1
  o = M;
else
  o = mean(M);
end
