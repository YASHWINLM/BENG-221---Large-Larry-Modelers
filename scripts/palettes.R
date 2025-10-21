library('ggplot2')
library('RColorBrewer')
library('khroma')
library('colorspace')

## Generate light and dark variants of country
## Tol's high contrast palette, 3-colors
# https://sronpersonalpages.nl/~pault/data/colourschemes.pdf

## Statistical significance
col_pal_signif = c('p < 0.05 *' = '#4477AA', 'p < 0.1 .' = '#66CCEE', 'N.S.' = '#BBBBBB',
                   'q < 0.05 *' = '#4477AA', 'q < 0.1 .' = '#66CCEE', 'N.S.' = '#BBBBBB')
col_pal_signif_pale = col_pal_signif %>%
  lapply(function(co) hex(mixcolor(alpha = 0.5, color1 = hex2RGB(co), color2 = hex2RGB('#FFFFFF')))) %>%
  do.call('c',.)

scale_fill_signif = scale_fill_manual(name = 'Significance', breaks = names(col_pal_signif), values = col_pal_signif)
scale_color_signif = scale_color_manual(name = 'Significance', breaks = names(col_pal_signif), values = col_pal_signif)
