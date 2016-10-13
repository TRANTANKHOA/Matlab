%% EXPAND PROCEDURE
% This is a private script for expanding the slice
%% Generate a random initial interval [a,b] for each sample
particles.b      = rand(target.M,1);
particles.a      = particles.b - 1;
% Marking the status of expansion for each sample. The expansion will
% continue for one sample if a flag is true (1) and stop if otherwise (0).
particles.flags  = true(target.M,1);    
%% Expanding the upper bound of the inteval [a,b]
m_w = 100;                  % Fix the maximum number of expansions
m_b = randi(m_w+1,target.M,1)-1;   % Allocate the maximum number of expansions of the upper bound
i=0;                        % Reset the counter for the number of expansions
while any(particles.flags)
    particles.pros(particles.flags,:)   = target.map(particles.pos(particles.flags,:),...
                                        particles.dirc(particles.flags,:),particles.b(particles.flags));  % x = theta+v.b
    particles                           = target.set(particles);
    particles.flags(particles.flags)    = particles.flags(particles.flags) & (i<m_b(particles.flags));      % is i<m_b ?
    particles.b(particles.flags)        = particles.b(particles.flags)+1;         % b = b+1
    i = i+1;
end

%% Expanding the lower bound of the inteval [a,b] (similar as above)
particles.flags = true(target.M,1);
m_a = m_w - m_b;            % Allocate the maximum number of expansions of the lower bound
i=0;    
while any(particles.flags)
    particles.pros(particles.flags,:)   = target.map(particles.pos(particles.flags,:),...
                                        particles.dirc(particles.flags,:),particles.a(particles.flags));  % x = theta+v.a
    particles                           = target.set(particles);
    particles.flags(particles.flags)    = particles.flags(particles.flags) & (i<m_a(particles.flags));      % is i<m_a ?
    particles.a(particles.flags)        = particles.a(particles.flags)-1;         % a = a-1
    i = i+1;
end
%% Counting the number of density evaluation per slice
particles.count    = floor(particles.b) - ceil(particles.a)+2;
