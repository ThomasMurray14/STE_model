% for i = 1 : 100
eta1 = 1;
eta0 = 0;
muhat_0_set = [0:0.05:1];
u_0_set  = [0:0.2:1];



% play with different values of alpha 
% when alpha 0.02, no perceptual uncertainty, 
% when alpha 1, a lot of perceptual uncertainty and the input is simply the expectation

al_temp = .5;%
figure;
hold on
xlabel('u_0_set'); % Label x-axis
ylabel('mu1'); % Label y-axis
title(['al = ', num2str(al_temp)]);
for s = 1 :length(muhat_0_set)
%close all;
muhat_0 = muhat_0_set(s);
mu1 = get_mus(eta1, eta0, u_0_set, muhat_0, al_temp);
plot(u_0_set, mu1)
waitforbuttonpress
end
hold off



figure;
hold on
xlabel('muhat_0_set'); % Label x-axis
ylabel('mu1'); % Label y-axis
title(['al = ', num2str(al_temp)]);
for s = 1 :length(u_0_set)
%close all;
u_0 = u_0_set(s);
mu1 = get_mus(eta1, eta0, u_0, muhat_0_set, al_temp);
plot(muhat_0_set, mu1)

waitforbuttonpress
end
hold off


function mu1 = get_mus(eta1, eta0, u_0_set, muhat_0_set, al_temp)
for j = 1 : length(u_0_set)
    for i = 1: length(muhat_0_set)
        u_0 = u_0_set(j);
        muhat_0 = muhat_0_set(i);
        und1 = exp(-(u_0 -eta1)^2/(2*al_temp));
        und0 = exp(-(u_0 -eta0)^2/(2*al_temp));
        %         mu(k,1) = muhat(k,1) *
        mu1(i,j) =  muhat_0*und1/(muhat_0 *und1 +(1 -muhat_0) *und0);
        
    end
end
end

