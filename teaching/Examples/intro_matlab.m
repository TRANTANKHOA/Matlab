v = [1 2 3 4]
I = eye(4)
help eye
M = rand(4,4)
v*M
cos(v)
log(v)
which kalman
s = zeros(size(v))';
for i = 1:4
    clc
    s = s +v(i).*M(:,i)
    pause(1/4)
end
s2 = sum(bsxfun(@times,M,v),2)
clc
u = rand(2,2)
n = randn(1e5,1);hist(n)
close all
r = -exprnd(1,1e5,1);hist(r,50)
edit experiment_9
newfunc('name')
help name
newfunc('kalman')
which kalman
name(particles)
my.name(particles)
clear
%% Bonus, try this
n = 1e8;
tic
for i=1:n
    rand;
end
toc
tic
rand(n,1);
toc
