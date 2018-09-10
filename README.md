# reweightfes
R package for Bonomi and Tiwary reweighting of metadynamics and other biased molecular simulations 

Done:
* `read.colvar(filename, cvs=2:3, bias=4)`
* `+.colvar`
* `print.colvar`
* `summary.colvar`
* `plot.colvar`
* `weightgibbs(cvfile, npoints=60, maxfe=100, temp=300, eunits="kJ/mol")`
* `reweightbonomi(cvfile, npoints=60, maxfe=100, temp=300, eunits="kJ/mol")`
* `reweighttiwary(cvfile, hillsfile, npoints=60, maxfe=100, nfes=100, temp=300, gamma=10, eunits="kJ/mol")`

Todo:
* data
* `feprofgibbs(minima, colvar)`
* `feprofbonomi(minima, colvar)`
* `feproftiwary(minima, colvar)`
* `imin` for `reweighttiwary`?
* if placed to metadynminer, add option fes and minima without hills


