# reweightfes
R package for on-the-fly and Tiwary reweighting of metadynamics and other biased molecular simulations 

Done:
* `read.colvar(filename, cvs=2:3, bias=4)`
* `+.colvar`
* `print.colvar`
* `summary.colvar`
* `plot.colvar`

Todo:
* `weightboltzmann(cvfile, npoints=60, maxfe=100, temp=300, eunits="kJ/mol")`
* `reweightbonomi(cvfile, npoints=60, maxfe=100, temp=300, eunits="kJ/mol")`
* `reweightiwary(cvfile, hillsfile, npoints=60, maxfe=100, nfes=100, temp=300, gamma=10, eunits="kJ/mol")`


