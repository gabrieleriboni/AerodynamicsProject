function I = simpcomp(xvect, yvect)

h = (xvect(end) - xvect(1)) / length(xvect);

I = (h/6) * (yvect(1) + yvect(end) + 2 * sum(yvect(3:2:(end-2))) + 4 * sum(yvect(2:2:(end-1))));