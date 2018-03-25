% norm_sig - Normalizing the signal by baseline substraction.%
%
% Usage: 
%     [normsig,base_line]=norm_sig(signal)
%
% Input parameters:
%     signal    : original signal which is given as a ROW vector  
%
% Output parameters:
%     normsig   : normalized signal  
%     base_line : base_line of the signal 

function [normsig,base_line]=norm_sig(signal)
normsig=zeros(size(signal));
base_line=zeros(size(signal));

for i=1:1:size(signal,1)
    
    dx = diff(smooth(signal(i,:)));
    if abs(dx(1)) < 0.005 %&& abs(dx(end)) < 0.005
        %we can be quite sure, that we are not somewhere within a T wave ->
        %linear nomalizing makes sense
        slope=(signal(i,end)-signal(i,1))/(length(signal)-1);
    %     N = 10;
    %     slope = (median(signal(i, end-N:end)) - median(signal(i,1:N))) / (length(signal)-1);
        x=1:1:size(signal,2);
        base_line(i,:)=signal(i,1)+slope.*(x-1);
        normsig(i,:)=signal(i,:)-base_line(i,:);
        %normsig(i,:)=normsig(i,:)/max(abs(normsig(i,:)));
    else
        %it is better to do nothing here
        normsig(i,:) = signal(i,:);        
    end
end

        
        