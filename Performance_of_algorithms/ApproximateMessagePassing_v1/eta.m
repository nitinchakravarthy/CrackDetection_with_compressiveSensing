function y = eta(x, threshold)
% ETA performs a soft thresholding on the input x.
tau = threshold;
y =  max(0,1-tau./max(abs(x),1e-10)).*x;