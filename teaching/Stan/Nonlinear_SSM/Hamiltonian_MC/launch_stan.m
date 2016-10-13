%% Building the Stan model
model_code = {
'data {'
'    int<lower=1> T; // number of schools '
'    vector[T] y;      // observed output data'
'}'
'parameters {'
'    real theta; '
'    real phi;'
'    vector[T] alpha;'
'    vector[T] y_tilde;'
'}'
'transformed parameters {'
'    real tau_2;'
'    tau_2 <- exp(theta);'
'}'
'model {'
'    theta ~ normal(0, 1);'
'    phi   ~ normal(0, 1);'
'    alpha[1] ~ normal(0, sqrt(tau_2) );'
'for (t in 2:T) {'
'    alpha[t] ~ normal(exp(-phi*alpha[t-1]^2)*alpha[t-1], sqrt(tau_2) );'
'               }'
'for (t in 1:T) {'
'    y[t] ~ normal(phi*exp(alpha[t]), 1);'
'    y_tilde[t] ~ normal(phi*exp(alpha[t]), 1);'
'               }'
'}'
};

data = struct('T',T,'y',y');

model = StanModel('model_code',model_code);

%% Compile the Stan model
model.compile();

%% Create the StanFit object
iter = 2e4; warmup = 1e3; chains = 5; M = iter*chains;
fit = model.sampling('data',data,'warmup',warmup,'iter',iter,'chains',chains,'verbose',false);
% fit = stan('model_code',model_code,'data',data,'warmup',1e3,'iter',1e5,'chains',6,'verbose',true);

% blocking the command line
fit.block();
%% Plotting and cleaning
stan_exit;
drawnow