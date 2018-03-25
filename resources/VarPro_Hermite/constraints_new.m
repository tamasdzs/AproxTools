function [trans_lb,trans_ub,dilat_lb,dilat_ub] = constraints_new(p,rule,fs,Rpos,N)
%CONSTRAINTS_NEW Summary of this function goes here
%   Detailed explanation goes here
    %% Setting constraints for the parameters of the QRS, T, P waves.
%     avgQRSwidth=0.110 * fs;   %Average number of samples for the QRS complex. 
    minQRSwidth=0.03 * fs;   %Average number of samples for the QRS complex. Carl: Might be better to set it to a smaller value (e.g. 0.03)
%     avgTwidth=0.210   * fs;   %Average number of samples for the T wave. 
    minTwidth=0.042 * fs;
%     avgPwidth=0.115   * fs;   %Average number of samples for the P wave. 
    minPwidth=0.041   * fs;

    maxQRSwidth=0.215 * fs;   %Maximum number of samples for the QRS complex. 
%     maxQRSwidth=0.07 * fs;   %Maximum number of samples for the QRS complex. Carl: Might be better to set it to a smaller value (e.g. 0.15)
    maxTwidth=0.400   * fs;   %Maximum number of samples for the T wave. 
%     maxTwidth = 0.15 * fs;
    maxPwidth=0.190   * fs;   %Maximum number of samples for the P wave. 
%     maxPwidth=0.08   * fs;
    
    %The ith translation parameter is restricted to the interval [trans_lb(i),trans_ub(i)].
    %i=1,2,3 correspond to the translation parameter of the QRS, T, P waves, respectively.
    min_RT = 0.0919 * fs;    %minimal distance between Rpeak and Tpeak
    max_RT = 0.441 * fs;     %maximal distance between Rpeak and Tpeak
    
    min_PR = 0.051 * fs;     %minimal distance between P peak and Rpeak
    max_PR = 0.307 * fs;     %maximal distance between P peak and Rpeak
    
    
    trans_lb=[Rpos-minQRSwidth/2, Rpos+min_RT, Rpos-max_PR]; %Lower bounds for the translation parameters.
    trans_ub=[Rpos+minQRSwidth/2, Rpos+max_RT, Rpos-min_PR]; %Upper bounds for the translation parameters.
    
    %The ith dilation parameter is restricted to the interval [dilat_lb(i),dilat_ub(i)].
    %i=1,2,3 correspond to the dilation parameter of the QRS, T, P waves, respectively.
    dilat_lb=[2*rule/(sqrt(p(1))*maxQRSwidth),2*rule/(sqrt(p(2))*maxTwidth),2*rule/(sqrt(p(3))*maxPwidth)]; %Lower bounds for the dilation parameters.
    dilat_ub=[2*rule/(sqrt(p(1))*minQRSwidth),2*rule/(sqrt(p(2))*minTwidth),2*rule/(sqrt(p(3))*minPwidth)]; %Upper bounds for the dilation parameters.    

end

