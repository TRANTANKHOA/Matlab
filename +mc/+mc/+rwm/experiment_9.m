%%
% propos = @(v) log(normpdf(v,mu_v,sigma_v));
% x = 0; px = target(x);
% v = 0; pv = propos(v);

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