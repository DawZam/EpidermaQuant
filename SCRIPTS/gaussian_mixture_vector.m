function [pp_est,mu_est,sig_est] = gaussian_mixture_vector(data, Y, KS,quick_stop, draw)
% FUNCTION: 
% Gaussian mixture decomposition of 1D signal by EM
% algorithm. Initial conditions for EM found by dynamic programming
% approach. No. of components found by BIC criterion.
% INPUT:
% data - signal values [1xn]
% KS - max no. of components
% draw - if draw results
% OUTPUT:
% pp_est - vector of esimated weights
% mu_est - vector of estimated component means
% sig_est - vector of estimated standard deviations of components

% AUTHORS:
% Michal Marczyk and Andrzej Polanski
% emails: {michal.marczyk,andrzej.polanski}@polsl.pl

if min(size(data)) ~= 1
    error('data must be 1D signal.')
end

bin_edge_sum = sum(Y);

N = length(data);
BIC = nan(KS,1);
logL = BIC;
D = nan(1,KS); 
alpha = cell(1,KS); mu = alpha; sigma = alpha;
    
%histogram of input data (for drawing and IC)
% breaks = linspace(min(data), max(data), (min(max(20, round(sqrt(N))), 100)));
% [y,x] = hist(data, breaks);
% [y,x] = hist(data,max(min(20,round(sqrt(N))),100));
y = Y;
x = data;

% decomposition for 1 component
[alpha{1},mu{1},sigma{1},logL(1)] = EM_iter(data,1,mean(data),std(data),Y, N);
%BIC(1) = -2*logL(1) + 2*log(bin_edge_sum);   

%decomposition for >2 components
stop = 1;
k = 2;
Nb = length(x);
aux_mx = dyn_pr_split_w_aux(x,y);    
while stop && k < KS
%     if draw; disp(['k = ' num2str(k)]); end
    
    %find initial condition using DP approach
    [~,opt_part] = dyn_pr_split_w(x,y,k-1,aux_mx);
    part_cl=[1 opt_part Nb+1]; 
    pp_ini = zeros(1,k); mu_ini = zeros(1,k); sig_ini = zeros(1,k);
    for kkps=1:k
        invec = x(part_cl(kkps):part_cl(kkps+1)-1);
        yinwec = y(part_cl(kkps):part_cl(kkps+1)-1);
        wwec = yinwec/(sum(yinwec));
        pp_ini(kkps) = sum(yinwec)/sum(y);
        mu_ini(kkps) = sum(invec.*wwec);
        sig_ini(kkps) = 0.5*(max(invec)-min(invec));
    end
    
    %perform decomposition
    [alpha{k},mu{k},sigma{k},logL(k)] = EM_iter(data,pp_ini,mu_ini,sig_ini,Y, N);

    %check convergence
    BIC(k) = -2*logL(k) + (3*k-1)*log(bin_edge_sum);   
    D(k) = -2*logL(k-1) + 2*logL(k);

    if quick_stop
     if 1 - chi2cdf(D(k),3) > 0.2 && k > 5 
         stop = 0;
     end
    end
    k = k+1;
end
[~,cmp_nb] = min(BIC); 
% if draw; figure; plot(1:KS,BIC);xlabel('No. of components'); ylabel('BIC');end

pp_est = alpha{cmp_nb};
mu_est = mu{cmp_nb};
sig_est = sigma{cmp_nb};     

if draw
 
    f_temp = zeros(1e5,cmp_nb); 
    x_temp = linspace(min(data),max(data),1e5)';
    step_x_temp = [min(diff(x_temp)); diff(x_temp)];
    for k=1:cmp_nb; f_temp(:,k) = pp_est(k)*normpdf(x_temp,mu_est(k),sig_est(k)).*step_x_temp; end 
    
    figure; hold on; box on
    bar(x,y,[min(data) max(data)],'hist');
    plot(x_temp',mean(diff(x))*sum(y)*sum(f_temp./repmat(step_x_temp,1,cmp_nb),2),'r.');
    for a=1:cmp_nb
        plot(x_temp',mean(diff(x))*sum(y)*f_temp(:,a)./step_x_temp,'g');
    end
    title(['Gaussian  mixture decomposition with ' num2str(cmp_nb) ' components'])
    xlabel('Variable'); ylabel('Counts');
end
    
end %end function

function [pp_est,mu_est,sig_est,logL] = EM_iter(x,alpha,mu,sig, Y, N)

bin_edge_sum = sum(Y);

x = x(:); Y=Y(:)'; alpha = alpha(:)'; mu = mu(:)'; sig = sig(:)';
[x,ind] = sort(x);
Y = Y(ind);
sig2 = sig.^2;
change = Inf;   
count = 1;
eps_change = 1e-7;
KS = length(alpha);
% SW = ((max(x)-min(x))/(5*KS))^2;     %minimum variance
SW = 1;
while change > eps_change && count < 5000
    old_alpha = alpha;
    old_sig2 = sig2;
    
    f = zeros(KS,N); sig = sqrt(sig2);
    for a=1:KS
        f(a,:) = norm_pdf(x,mu(a),sig(a));
    end
    px = alpha * f;
    px(isnan(px) | px==0) = 5e-324;
    for a=1:KS
       pk = ((alpha(a)*f(a,:)).*Y)./px;
       denom = sum(pk);
       mu(a) = sum((pk*x)/denom);
       sig2num = sum(pk*((x-mu(a)).^2));
       sig2(a) = max(SW,sig2num/denom);
       alpha(a) = denom/bin_edge_sum;
    end
   
    change = sum(abs(alpha-old_alpha)) + sum(((abs(sig2-old_sig2))./sig2))/(length(alpha));
    count = count+1;
end

% RETURN RESULTS
logL = sum(log(px).*Y);
[mu_est,ind] = sort(mu);
sig_est = sqrt(sig2(ind));
pp_est = alpha(ind);

end %end function

function y = norm_pdf(x,mu,sigma)
y = exp(-0.5*(((x - mu)./sigma).^2)) ./(2.506628274631 * sigma);
end %end function

function aux_mx=dyn_pr_split_w_aux(data,ygreki)
N=length(data);
% aux_mx
aux_mx=zeros(N,N);
for kk=1:N-1
   for jj=kk+1:N
       aux_mx(kk,jj)= my_qu_ix_w(data(kk:jj-1),ygreki(kk:jj-1));
   end
end
end %end function

function [Q,opt_part]=dyn_pr_split_w(data,ygreki,K_gr,aux_mx)
% initialize
Q=zeros(1,K_gr);
N=length(data);
p_opt_idx=zeros(1,N);
p_aux=zeros(1,N);
opt_pals=zeros(K_gr,N);
for kk=1:N
    p_opt_idx(kk)=my_qu_ix_w(data(kk:N),ygreki(kk:N));
end
% iterate
for kster=1:K_gr
   for kk=1:N-kster
       for jj=kk+1:N-kster+1
           p_aux(jj)= aux_mx(kk,jj)+p_opt_idx(jj);
       end
       [mm,ix]=min(p_aux(kk+1:N-kster+1));
       p_opt_idx(kk)=mm;
       opt_pals(kster,kk)=kk+ix(1);
   end
   Q(kster)=p_opt_idx(1);
end
% restore optimal decisions
opt_part=zeros(1,K_gr);
opt_part(1)=opt_pals(K_gr,1);

if K_gr ~= 1
    for kster=K_gr-1:-1:1
       opt_part(K_gr-kster+1)=opt_pals(kster,opt_part(K_gr-kster));
    end
end
end %end function

function wyn=my_qu_ix_w(invec,yinwec)
PAR=1;
PAR_sig_min=0.1;
invec=invec(:);
yinwec=yinwec(:);
if (invec(length(invec))-invec(1))<=PAR_sig_min || sum(yinwec)<=1.0e-3
    wyn=inf;
else
   wwec=yinwec/(sum(yinwec));
   wyn1=(PAR+sqrt(sum(((invec-sum(invec.*wwec)).^2).*wwec)))/(max(invec)-min(invec));
   wyn=wyn1;
end
end %end function