%% SHRINK PROCEDURE
% This is a private script for sampling from the slice
particles.flags = true(target.M,1);
rnd = zeros(target.M,1);
while any(particles.flags)
    particles.count(particles.flags)    = particles.count(particles.flags) + 1;
    rnd(particles.flags)                = (particles.b(particles.flags)-particles.a(particles.flags)).*...
                                        rand(nnz(particles.flags),1)+particles.a(particles.flags);
    particles.pros(particles.flags,:)   = target.map(particles.pos(particles.flags,:),...
                                        particles.dirc(particles.flags,:),rnd(particles.flags));  % x = theta+v.r 
    flags                               = particles.flags;
    particles                           = target.set(particles);                         
    particles.flags(flags)              = ~particles.flags(flags);
    idx                                 = particles.flags;
    idx(particles.flags)                = rnd(particles.flags)>0;
    particles.b(idx)                    = rnd(idx);
    particles.a(~idx & particles.flags) = rnd(~idx & particles.flags);
end
particles.pos  = particles.pros;
