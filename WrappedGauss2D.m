% y=WrappedGauss2D(x,y,mu,sigma,k)
% x and y are the corresponding x and y arguments to the 2D Gaussian
% and must be the same size.  mu and sigma are the mean and std, must be
% numbers. k is the number of times to wrap.  in most cases, k=2 or 3 is
% fine.

function g=WrappedGauss2D(x,y,mu,sigma,k)

    [m1,m2]=size(x);
    [n1,n2]=size(y);
    if(m1~=n1 || m2~=n2 || numel(mu)~=1 || numel(sigma)~=1 || numel(k)~=1)
        error('x and y must be same size; mu and sigma must be numbers.')
    end

    g=zeros(size(x));
    for j1=-k:k
    for j2=-k:k    
        g=g+normpdf(x+j1,mu,sigma).*normpdf(y+j2,mu,sigma);
    end
    end

end

