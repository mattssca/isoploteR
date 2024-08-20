<!-- badges: start -->
<!-- badges: end -->

# isoploteR
Wrangles expression matrix into shape, calculate isoform frquencies and ploting.

## Installation
Run the following command to install the development version from GitHub:
```
devtools::install_github("mattssca/isoploteR")
```

## Example
``` r
my_samples = colnames(expression_sub)[-1]
my_samples = my_samples[-1]


my_samples = colnames(expression_sub)[-1]

#run function
these_isos = get_isos(this_data = expression_sub,
                      annotations = gene_annotations,
                      these_sample_ids = my_samples,
                      these_isoforms = c("ENST00000395080",
                                         "ENST00000237623",
                                         "ENST00000360804",
                                         "ENST00000508233",
                                         "ENST00000681973"),
                      plot_title = "SPP1",
                      plot_subtitle = "Isoforms Frequency")
                      
```
<img src="/figs/example_plot.png" width="100%" />
