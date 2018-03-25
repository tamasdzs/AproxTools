%% Testing the correctness of the derivated Hermite functions.
basenums=4;
x=linspace(-5,5,300);
[H, DH] = derivated_hermite_system(x,basenums);
%Computing divided differences with stepsize h=x(2)-x(1).
divdiff=diff(H)/(x(2)-x(1));
figure(1);
plot(DH,'r'); hold on; plot(divdiff,'b'); hold off;
title('Analytic derivatives (red), divided differences (blue)');

%% Testing the correctness of the partial derivatives with respect to the free parameters: dilations, translations.
y=1:300; N=300;
ada=@(alpha) ada_Hermite(y,basenums,alpha,false);
alpha=[0.1, 10];
[Phi,dPhi]=ada(alpha);
h=[0.001,0];
[Phi_h,dPhi_h]=ada(alpha+h);

figure(2);
subplot(1,2,1);
divdiff=(Phi_h-Phi)/h(1);
plot(dPhi(:,1:2:end),'r'); hold on; plot(divdiff,'b'); hold off;
title('For dilation: Analytic derivatives (red), divided differences (blue)');

subplot(1,2,2);
h=[0,0.1];
[Phi_h,dPhi_h]=ada(alpha+h);
divdiff=(Phi_h-Phi)/h(2);
plot(dPhi(:,2:2:end),'r'); hold on; plot(divdiff,'b'); hold off;
title('For translation: Analytic derivatives (red), divided differences (blue)');
