%% Buckling Mode Prediction (contour)
clear
close all

nu = 0.5; % Possion's ratio of rubber
alpha = 0.5:-0.001:0.1; % Radius/Height
beta = 0.01:0.0001:0.05; % thickness/Height
n=[2,3,4,5,6,7,8]; % initial guess of fold numbers

p_cr = zeros(length(alpha), length(beta));
N = zeros(length(alpha), length(beta));
for p = 1:length(alpha)
    aa = alpha(p);
    for q = 1:length(beta)
        bb = beta(q);
        p_n = zeros(1,length(n));
        for i = 1:length(n)
            p_n(i) = ((pi*aa)^2+n(i)^2)^2*(bb/aa)^3/(12*(1-nu^2)*n(i)^2)+...
                (pi*aa)^4*(bb/aa)/(n(i)^2*((pi*aa)^2+n(i)^2)^2);
                %Z(i)=(1/aa*1/bb)*sqrt(1-nu^2);
            %p_n(i) = (1+n(i)^2)^2/n(i)^2+(12/pi^4)*Z(i)^2/(n(i)^2*(1+n(i)^2)^2);
        end
        [p_cr(p,q),I] = min(p_n);
        N(p,q) = n(I);
    end
end
%M = flipud(N)
figure
%h = heatmap(beta,alpha,N,'CellLabelColor','none')
%contourf(beta,alpha,N,'--');
h = surf(beta,alpha,N,'FaceAlpha',1);
h.EdgeColor = 'none';
set(gca,'FontSize',16)
pbaspect([1 1 1])
%grid off
shading interp;
colormap(parula(8));
view(0,90)
xlabel('t/H','FontSize',16);
ylabel('R/H','FontSize',16);

%find_n_alt(0.2, 0.025)
n_calc(0.025, 0.2)

function n_find = n_calc(aa,bb)
n=[2,3,4,5,6,7,8]; % initial guess of fold numbers
nu = 0.5;
p_n = zeros(length(n),1);
for i = 1:length(n)
    p_n(i) = ((pi*aa)^2+n(i)^2)^2*(bb/aa)^3/(12*(1-nu^2)*n(i)^2)+...
        (pi*aa)^4*(bb/aa)/(n(i)^2*((pi*aa)^2+n(i)^2)^2);
        %Z(i)=(1/aa*1/bb)*sqrt(1-nu^2);
    %p_n(i) = (1+n(i)^2)^2/n(i)^2+(12/pi^4)*Z(i)^2/(n(i)^2*(1+n(i)^2)^2);
end
[~,I] = min(p_n);
n_find = n(I);
end