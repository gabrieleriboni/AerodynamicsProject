function I = trapcomp(xvect, yvect)

h = (xvect(end) - xvect(1)) / length(xvect);

I = h/2 * (2 * sum(yvect) - yvect(1) - yvect(end));