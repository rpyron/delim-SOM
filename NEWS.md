# delimSOM 2.1.0

## Changed

* **Breaking change.** Training, clustering, and leave-one-layer-out analysis
  now load saved results only on explicit request. Existing `.Rdata` files
  remain readable. To reuse one, give its path and set `reuse.SOM.results = TRUE`
  for training or clustering. For leave-one-layer-out analysis, set
  `reuse.leave.one.layer.out.results = TRUE` and
  `save.leave.one.layer.out.results = FALSE`. To save fresh results at an
  existing path, explicitly allow overwrite or choose a new path.
* Pooled species coefficients now align specimen-group labels across occupied
  K before counting assignments. Each retained map still contributes one vote.
  Coefficients can change when maps have different numbers of occupied groups.
* Specimen IDs must be unique, nonblank, and non-missing within each layer.
  Unnamed single-layer inputs and compatible older single-layer fits can still
  use positional IDs.

## Fixed

* Match `kohonen`'s specimen-to-unit distances for soft assignments and
  topographic error, including missing-value scaling and the Tanimoto boundary.
* Preserve soft weights for unoccupied grid clusters and retain valid maps
  whose specimen groups can be relabeled even when the grid contains more
  clusters.
* Calculate assignment diagnostics before padding. Return `NA` when a
  one-cluster map makes a diagnostic undefined, including leave-one-layer-out
  margin comparisons.
* Keep automatic grids within the retained specimen count. Estimate single-layer
  grid shape from the two leading axes, including when there are more loci than
  specimens.
* Render unoccupied grid clusters in gray and use stored specimen labels for
  occupied-cluster colors, including when plotting older saved fits.

## Added

* Expanded cross-platform functional and regression tests and included the
  GPL-3 license text.

# delimSOM 2.0

* Converted delim-SOM 2.0 into an installable R package.
* Added package documentation and exported public functions.
* Added initial package API tests.
