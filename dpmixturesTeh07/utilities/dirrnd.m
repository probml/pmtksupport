function pi = dirrnd(aa);
% pi = dirrnd(aa)
% draws a sample from a dirichlet with parameter vector aa

pi = randg(aa);
pi = pi/sum(pi);
