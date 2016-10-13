%% Building the Stan model
model_code = {
'data {'
'    int<lower=1> T; // number of schools '
'    real y[T];      // observed output data'
'}'
'parameters {'
'    real theta; '
'    real y_tilde[T];'
'}'
'transformed parameters {'
'    real tau_2;'
'    tau_2 <- exp(theta);'
'}'
'model {'
'    theta ~ normal(0, 1);'
'    y ~ normal(0, sqrt(1 + tau_2) );'
'    y_tilde ~ normal(0, sqrt(1 + tau_2) );'
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