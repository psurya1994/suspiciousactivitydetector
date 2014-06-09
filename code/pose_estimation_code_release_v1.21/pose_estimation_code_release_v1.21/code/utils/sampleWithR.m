function Y = sampleWithR(weights,K)
%Y = sampleWithR(weights,K)
%Generates K samples from the discrete distribution specified by weights
%O(K) 

%Need to handle border effects
cdf = cumsum(weights);
Y = histc(rand(K,1)*cdf(end),[0; cdf]);
Y = [Y(1:end-2); Y(end-1)+Y(end)];
