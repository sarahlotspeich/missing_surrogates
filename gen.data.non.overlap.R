# Functions to generate data 
gen.data.non.overlap = function(setting, n1, n0) {
  s1 = g.1(n1)
  y1 = f.cond.1(s1)
  s0 = g.0(n0)
  y0 = f.cond.0(s0)
  return(data.frame("s1" = s1, "y1" = y1, "s0" = s0, "y0" = y0))
}
f.cond.1 = function(s.vector) {
  eps1 = rnorm(length(s.vector),0,3)
  y1 = 2+5*s.vector+1 + 1*s.vector + eps1
  return(y1)		
}
f.cond.0 = function(s.vector) {
  eps0 = rnorm(length(s.vector),0,3)
  y0 = 2+5*s.vector+ eps0
  return(y0)		
}
g.1 = function(n, alpha0=5) { return(rnorm(n, alpha0 + 1,1/2))}
g.0 = function(n,alpha0=5) { return(rnorm(n, alpha0,1))}