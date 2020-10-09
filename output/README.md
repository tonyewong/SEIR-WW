# Output directory

The `infections_failurecaseX_25Sep2020.csv` files give the number of overall cumulative infections for the
entire semester under the different failure cases and different combinations of screening test parameters
$T_{s,ww}$ and $T_{lag}$. The `25Sep2020` files consider $T_{lag}$ values of 1-6 days, and $T_{s,ww}$ values of
1, 2, 4, 6, and 8 days. The files from `05Oct2020` are the same, except $T_{lag}$ is held at 1 day and
$T_{s,ww}$ values of all integers 1-8 days are considered.

The file `risk_results_25Sep2020.RData` includes the $T_{lag}$ parameter at values of 1, 2, 3, 4, 5, and 6
days of delay before initiating the wastewater-triggered screening tests, and $T_{s,ww}$ at values of 1, 2, 4,
6, and 8 days to screen the entire wastewater-contributing subpopulation. The `05Oct2020` RData results file
corresponds to the $T_{lag}=1$ day version.

The three `sobol_...rds` files correspond to the results from the Sobol' sensitivity analyses. The designation
of `priors` within the file names (`high`, `mid` (control), or `small`) indicates whether the file corresponds
to the large university case, RIT control case, or small-college case, respectively.
