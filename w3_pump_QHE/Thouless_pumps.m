%-------------------1D Thouless pumps---------------
%-----H = k^2/2m + A(1-cos(x/lambda))
clc;
clear;
close all
%----------------------------------------
n = 15;     %the number of time iterate wire
L = n*17;
k =-pi:2*pi/L : (L-2)*pi/L;
a = 1; 
omega = 2*n*pi/L;
mu = 0;
t = 1;
A = 0.8;
%-----------------------------
% e_A = zeros(size(k,2),size(k,2));
%------------------------------
U = @(i) repmat(exp(1i*k(i)*linspace(1,L,L)*a).',1,L);
H = diag((2*t-mu + 1*A)*ones(1,L)) - diag(A*cos(omega*linspace(1,L,L))) - diag(t*ones(1,L-1),-1)- diag(t*ones(1,L-1),1);
% H = diag((2*t-mu )*ones(1,L)) - diag(A*cos(omega*linspace(1,L,L))) - diag(t*ones(1,L-1),-1)- diag(t*ones(1,L-1),1);
H(1,L) = -t;
H(L,1) = -t;
H_k = round(fft2(H),10)./L;                   %2-D fourier transfrom
H_k = fftshift(H_k);                          %shift zero frequency component to center of spectrum
H_k = flip(H_k);                              %shift to block diag
[V,e]=eig(H_k);
e_A = diag(e);
e_A = [flip(e_A(1:2:end));e_A(2:2:end)];      %separate odd and even indice 

E=reshape(e_A,n,[]);                          %reshape to correct brillouin zone
E = [flip(E(1:round(n/2),:),2);E(round(n/2)+1:end,:)];
figure(1);hold on;
k_wire = -pi:2*pi/(n) : pi-2*pi/(n);
for ii = 1:size(E,2)
plot(k_wire,E(:,ii))
end
xlim([-pi,pi])
ylim([-0.2,1.3])
xlabel('ka')
ylabel('engergy [t]')
figure(2)
plot(k,e_A)
% line([ -pi,pi],[0,0],'linestyle','--','color','r');
