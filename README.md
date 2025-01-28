General Description:
This analysis begins with the script START_run_tortuosity_routines_make_results.m This will load pre-processed data structs and then go through the analyses. david_sinuosity3.m is used to find behavioural regions where animals move from "search" to "exploit" and vice versa. This will also put out "control" regions where no transitions occur. The script steps through some analyses for low-to-high tortuosity segments and high-to-low tortuosity segments. These results are stored in structs, which are then saved all together. The analyses are:
Calculation of event-related potential and phase,
Spectral coherence,
Phase lag index,
Power correlations,
Cross-spectrum coupling.
The resultant saved structs are run by analyse_and_plot_tortuosity_results.m which will do some data munging and then plot/analyse the results.
