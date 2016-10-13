%% Building the Stan model
model_code = {
'functions {'
'// ... function declarations and definitions ...'
'}'
'data {'
'    int<lower=1> T; // number of data '
'    real y[T];      // observed output data'
'}'
'parameters {'
'    real phi;'
'    real logsigma; '
'    real logtau;'
'    real x[T];'
'    real y_tilde[T];'
'}'
'transformed parameters {'
'    real tau;'
'    real sigma;'
'    real kappa;'
'    tau   <- sqrt(exp(logtau));'
'    sigma <- sqrt(exp(logsigma));'
'    kappa <- tau/sqrt(1-phi^2);'
'}'
'model {'
'    logtau   ~ normal(0, 1);'
'    logsigma ~ normal(0, 1);'
'    phi   ~ normal(0, 1);'
'    x[1]  ~ normal(0,kappa);'
'for (t in 2:T) {'
'    x[t] ~ normal(phi*x[t-1], tau);'
'               }'
'    y       ~ normal(x, sigma);'
'    y_tilde ~ normal(x, sigma);'
'}'
};

data = struct('T',sys.T,'y',sys.y');

model = StanModel('model_code',model_code);

%% Compile the Stan model
model.compile();

%% Create the StanFit object
iter = 2e3; warmup = 1e3; chains = 5; M = iter*chains;
fit = model.sampling('data',data,'warmup',warmup,'iter',iter,'chains',chains,'verbose',false);
% fit = stan('model_code',model_code,'data',data,'warmup',1e3,'iter',1e5,'chains',6,'verbose',true);

% blocking the command line
fit.block();
%% Plotting and cleaning
stan_exit;
drawnow