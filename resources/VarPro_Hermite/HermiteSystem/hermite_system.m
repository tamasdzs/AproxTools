%% Computing the matrix of Hermite functions (N x n). 
function H = hermite_system(x,n)
    x=reshape(x,length(x),1);
    H=hermite_pol(n,x);
    w=exp(-x.^2/2);                  %weight function;
    for l=1:n
        H(:,l)=w.*H(:,l)./sqrt(2^(l-1)*factorial(l-1)*sqrt(pi));
    end
    H=H(:,1:n);
end

function H  = hermite_pol( n, x )
    if 0==n
        H=ones(length(x),1);         %zero order Hermite polynomial
    else    
        H=zeros(length(x),n);
        %weight function
        H(:,1)=1;                    %zero order Hermite polynomial
        H(:,2)=2*x;                  %first order Hermite polynomial

        %Hermite polynomials by recursion.
        for i=2:1:n-1
            H(:,i+1)=2*(x.*H(:,i)-(i-1)*H(:,i-1));
        end
    end
end