function A = local_sum_zero(A,m,n)
%A = local_sum_rep(A,m,n)
%Computes all the m X n sums for the local_sum
%Zero pads by replicating the border of A
%We'll only return the portion which is valid
%[mm,nn] =size(A);

%B = [zeros(m,nA+2*n-1,dZ);
%     zeros(mA,n,dZ) A zeros(mA,n-1,dZ);
%     zeros(m-1,nA+2*n-1,dZ)];
%B = [zeros(m,nA+n,dZ);zeros(mA,n,dZ) A];
%s = cumsum(B,1);
[mm,nn] = size(A(:,:,1));
if nargin == 2,
    n = m;
end
if m == 0 && n == 0
  return
end
A = padarray(A,round([m n]/2),'both');
s = cumsum(A,1);
c = s(1+m:end,:,:)-s(1:end-m,:,:);
s = cumsum(c,2);
A = s(:,1+n:end,:)-s(:,1:end-n,:);
A = A(1:mm,1:nn,:);
