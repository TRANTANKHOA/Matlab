%% Building the Stan model
model_code_full = {
'data {'
'    int<lower=1> T; // number of schools '
'    real y[T];      // observed output data'
'}'
'parameters {'
'    real theta; '
'    real alpha[T];'
'    real y_tilde[T];'
'}'
'transformed parameters {'
'    real tau_2;'
'    tau_2 <- exp(theta);'
'}'
'model {'
'    alpha ~ normal(0, 1);'
'    theta ~ normal(0, 1);'
'    y ~ normal(alpha, sqrt(tau_2) );'
'    y_tilde ~ normal(alpha, sqrt(tau_2) );'
'}'
};

data = struct('T',T,'y',y');

model_full = StanModel('model_code',model_code_full);

%% Compile the Stan model
model_full.compile();

%% Create the StanFit object
fit_full = model_full.sampling('data',data,'warmup',warmup,'iter',iter,'chains',chains,'verbose',false);
% fit_full = stan('model_code',model_code,'data',data,'warmup',1e3,'iter',1e5,'chains',6,'verbose',true);

% blocking the command line
fit_full.block();
