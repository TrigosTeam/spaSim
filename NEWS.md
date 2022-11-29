# spaSim 1.1.2

BUG FIXES

* Fixed documentation of argument `jitter` in `TIS()` and `simulate_background_cells()`.
* Fixed incorrectly plotted image in README.

# spaSim 1.1.1

SIGIFICANT USER CHANGE

* `simulate_background_cells()` added an option (`method = "Even"`) to simulate
evenly spaced background images. This change is accompanied with addition of two 
parameters, `method` - to choose the background cell distribution and `jitter`
- the parameter to simulate evenly spaced background cells. Tutorials and 
function `TIS()` were modified accordingly.

# spaSim 0.99.5

SIGNIFICANT USER CHANGE

* `TIS()` only plots one plot when arg `plot_image` is `TRUE`.

# spaSim 0.99.4

SIGNIFICANT USER CHANGE

* Added an optional parameter `plot_image` to `simulate_background_cells()`.

# spaSim 0.99.3

NO SIGNIFICANT CHANGE

# spaSim 0.99.2

NEW FEATURES

* Added `utils.R` file to store internal functions.

SIGNIFICANT USER-VISIBLE CHANGES

* Updated the main object class from `SingleCellExperiment` to `SpatialExperiment`.
* Changed the parameter `Phenotype` to `Cell.Type`.

# spaSim 0.99.1

NEW FEATURES

* Added parameter validity checks.

SIGNIFICANT USER-VISIBLE CHANGES

* Improved parameter descriptions.

BUG FIXES

* None

# spaSim 0.99.0

NEW FEATURES

* This is the first version of spaSim.
