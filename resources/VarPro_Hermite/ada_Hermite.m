%% This variation of the VarPro algorithm was made for least square fitting of dilated and translated 
%% Hermite functions. In [2], we adopted the original algorithm of [1] to these Hermite functions. 
%% Namely, the values of the Hermite functions are stored in matrices and their 
%% pseudoinverses are calculated via economic singular value decomposition.   

%% [1] Dianne P. O'Leary and Bert W. Rust, Variable Projection for Nonlinear Least Squares Problems,
%%     Computational Optimization and Applications (2012), doi 10.1007/s10589-012-9492-9.
%% [2] ???.

%  Last Modified: March 15, 2018.
%  Version 1.0.
%
% ada_Hermite - Implementing the 'ada' function for Hermite functions with dilation and translation parameters.
%               This function is for general use in least square fitting of Hermite functions with free parameters.           
%
% Usage: 
%     [Phi,dPhi,Ind] = ada_Hermite(y,alpha,show);
%
% Input parameters:
%     y       : ordinates of the sample points. 
%     alpha   : vector of free parameters, i.e. dilation and translation parameter of the Hermite functions.
%     show    : display the least fitting of Hermite functions at each step of the algorithm (true/false).
%
% Output parameters:
%     Phi  : The ith column contains the samples of the ith basis functions. 
%     dPhi : The columns contain the partial derivative information for Phi. 
%     Ind  : Column l of dPhi contains the partial derivative of Phi(:,j) with respect to alpha(i), 
%            where j = Ind(1,l) and i = Ind(2,l). Those partial derivatives that are independent of 
%            alpha(i) are equal to zero and so they are not stored in dPhi. That is why the index array 'Ind' is required.
%
%  Copyright (c) 2018, Péter Kovács <peter.kovacs@jku.at> (kovika@inf.elte.hu)  
%  Johannes Kepler University, Linz, Austria, 2018.   
%   
%  Permission to use, copy, modify, and/or distribute this software for  
%  any purpose with or without fee is hereby granted, provided that the  
%  above copyright notice and this permission notice appear in all copies.  
%  
%  THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL  
%  WARRANTIES WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED  
%  WARRANTIES OF MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR  
%  BE LIABLE FOR ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES  
%  OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS,  
%  WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION,  
%  ARISING OUT OF OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS  
%  SOFTWARE. 

function [Phi,dPhi,Ind] = ada_Hermite(y,hermitebasenum,alpha,show)
%
%     This function computes Phi and, if possible, dPhi.
%
%     On Input: 
%
%        alpha q x 1    contains the current value of the alpha parameters.
%
%        Note:  If more input arguments are needed, use the standard
%               Matlab syntax to accomplish this.  For example, if
%               the input arguments to ada are t, z, and alpha, then
%               before calling varpro, initialize t and z, and in calling 
%               varpro, replace "@ada" by "@(alpha)ada(t,z,alpha)".
%
%     On Output:
%
%        Phi   m x n1   where Phi(i,j) = phi_j(alpha,t_i).
%                       (n1 = n if there is no extra term; 
%                        n1 = n+1 if an extra term is used)
%        dPhi  m x p    where the columns contain partial derivative
%                       information for Phi and p is the number of 
%                       columns in Ind 
%                       (or dPhi = [] if derivatives are not available).
%        Ind   2 x p    Column k of dPhi contains the partial
%                       derivative of Phi_j with respect to alpha_i, 
%                       evaluated at the current value of alpha, 
%                       where j = Ind(1,k) and i = Ind(2,k).
%                       Columns of dPhi that are always zero, independent
%                       of alpha, need not be stored. 
%        Example:  if  phi_1 is a function of alpha_2 and alpha_3, 
%                  and phi_2 is a function of alpha_1 and alpha_2, then 
%                  we can set
%                          Ind = [ 1 1 2 2
%                                  2 3 1 2 ]
%                  In this case, the p=4 columns of dPhi contain
%                          d phi_1 / d alpha_2,
%                          d phi_1 / d alpha_3,
%                          d phi_2 / d alpha_1,
%                          d phi_2 / d alpha_2,
%                  evaluated at each t_i.
%                  There are no restrictions on how the columns of
%                  dPhi are ordered, as long as Ind correctly specifies
%                  the ordering.
%
%        If derivatives dPhi are not available, then set dPhi = Ind = [].
N=length(y);
if ~mod(N,2)
    t=(-N/2:N/2-1)';
else
    t=(-floor(N/2):floor(N/2))';
end
dilation=alpha(1); translation=round(N/2)-alpha(2);

x=dilation*(t+translation);
[H, DH] = derivated_hermite_system(x,hermitebasenum);
Phi=H;
dPhi=zeros(N,2*hermitebasenum);
Ind=zeros(2,2*hermitebasenum);
Ind(1,1:2:end)=1:hermitebasenum;
Ind(1,2:2:end)=1:hermitebasenum;

dPhi(:,1:2:end)=DH.*(t+translation); %Derivatives with respect to the dilation parameter.
Ind(2,1:2:end)=ones(1,hermitebasenum);

dPhi(:,2:2:end)=-DH*dilation; %Derivatives with respect to the translation parameter.
Ind(2,2:2:end)=2*ones(1,hermitebasenum);

%Displaying the results at each step.
if show
    aprx=Phi*(Phi\y); 
    prd=norm(y-reshape(aprx,size(y)))/norm(y-mean(y))*100;
    display(sprintf('PRD of the VarPro approximation: %.2f%%',prd));
    plot(1:N,y,'b',1:N,aprx,'r','LineWidth',2);
    legend('ECG',sprintf('VarPro PRD: %.2f%%',prd));
    drawnow;
    pause(0.5);
end