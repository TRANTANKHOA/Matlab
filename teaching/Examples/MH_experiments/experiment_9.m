clear; clc; n=1e5; particles = zeros(n,2); avg_acc = 0;
%%
mu_x = 1; sigma_x = 1.5;
target = @(x) log(normpdf(x,mu_x,sigma_x));
mu_v = 0; sigma_v = 1.55*sqrt(2.38)*sigma_x;
%%
propos = @(v) log(normpdf(v,mu_v,sigma_v));
x = 0; px = target(x);
v = 0; pv = propos(v);

for i=1:n
    %% update h
    log_h = px + pv + log(rand);
    %% update v
    rho = sqrt(-2*(log_h - px) - log(2*pi*sigma_v^2))*sigma_v;
    v = (rand*2-1)*rho + mu_v; pv = propos(v);
    %% update v and x
    v_new = -v;     pv_new = propos(v_new);
    x_new = x + v;  px_new = target(x_new);
    %% Accept - reject test
    if (px_new + pv_new) >= log_h
        x = x_new; px = px_new;
        v = v_new; pv = pv_new;
        acc = 1;
    else
        acc = 0;
    end
    avg_acc = avg_acc*(i-1)/i + acc/i;
    particles(i,:) = [x,v];
end
keyboard
plot_time_series(particles(:,1),'x',20,'Slice Kernel')
subplot(2,2,2)
hold on
y = -3*sigma_x:0.1:3*sigma_x;
plot(normpdf(y,mu_x,sigma_x),y)