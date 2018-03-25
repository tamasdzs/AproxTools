%Demonstrating varpro for Bspline ECG knot optimization

%Loading data, libraries, etc.
addpath('./HermiteSystem');

% Example for the optimization process.
load('beat1'); Rpos=117; cut=125; amplitude=0;
%load('beat2'); Rpos=136; cut=155; amplitude=0.0;
N=length(beat);
baseline=(1:N)';
baselineshift=amplitude*[ones(cut,1);-ones(N-cut,1)];
signal=reshape(beat,N,1)+baselineshift; %signal should be a column vector.
orig_sig=signal;

% Setting parameters.
show=true;
options = optimset('lsqnonlin');
options = optimset(options,'MaxIter',50);
w=ones(size(signal)); %Uniform weighting.

basenums=[4 4 4]; 
numfunsys=2; %If you extend the number of function systems used in the expansion, set the proper number of free parameters as well.
x=zeros(3,numfunsys);
co=cell(3,numfunsys);

% Making constraints for the heartbeat signal.
p=[2 2 2]; rule=2; fs=250;
[trans_lb,trans_ub,dilat_lb,dilat_ub] = constraints_new(p,rule,fs,Rpos,N);
%% Optimization
for i=1:1:length(basenums)
    %% Constraints for Hermite functions
    lb=[dilat_lb(i);trans_lb(i)];   %Constraints for the parameters of the Hermite functions. (LB)
%    lb=[lb; 1e-1; trans_lb(i)];     %Constraints for the parameters of the Sigmoid function.  (LB)
    ub=[dilat_ub(i);trans_ub(i)];   %Constraints for the parameters of the Hermite functions. (UB)
%    ub=[ub; 1e+2; trans_ub(i)];     %Constraints for the parameters of the Sigmoid function.  (UB)

    %ada function for evaluating basis functions and their partial derivatives.
    ada=@(alpha) ada_Hermite(signal,basenums(i),alpha,show);
    n=basenums(i);
    
    %even coordinates represent the translation parameters
    %odd coordinates represent the dilation parameters
    x0=(lb+ub)/2;

    [x(i,:), co{i}, wresid, wresid_norm, aprx, Regression] = varpro(signal, w, x0, n, ada, lb, ub, options);
    %x(i,:)=alpha; x(i,2:2:end)=round(N/2)-x(i,2:2:end);
    signal=signal-aprx;    
end

%% Reconstruction 
aprx=zeros(length(basenums),N);
for i=1:1:length(basenums)
    Phi=ada_Hermite(signal,basenums(i),x(i,:),false);
    aprx(i,:)=Phi*co{i};
end

%% Displaying the results
subplot(2,1,1);
plot(1:N,orig_sig,'b',1:N,sum(aprx),'r');
legend('Original ECG','Approximated');
subplot(2,1,2);
plot(1:N,aprx(1,:),'r',1:N,aprx(2,:),'g',1:N,aprx(3,:),'b');
legend('QRS','T','P');


%% Displaying the approximation error
PRD=norm(orig_sig-sum(aprx)')/norm(orig_sig-mean(orig_sig))*100;
display(sprintf('PRD: %.2f%%',PRD));
