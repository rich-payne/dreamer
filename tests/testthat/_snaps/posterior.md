# continuous predictive runs

    Code
      stats
    Output
      # A tibble: 3 x 5
         dose reference_dose  mean `2.50%` `97.50%`
        <dbl>          <dbl> <dbl>   <dbl>    <dbl>
      1     1              0  1.99  -0.344     3.98
      2     3              0  5.54   4.17      6.82
      3     5              0  9.43   7.62     11.6 

# binary predictive runs

    Code
      stats
    Output
      # A tibble: 3 x 5
         dose reference_dose  mean `2.50%` `97.50%`
        <dbl>          <dbl> <dbl>   <dbl>    <dbl>
      1     1              0  0.16  0         0.478
      2     3              0  0.32  0.0225    0.655
      3     5              0  0.22  0         0.555

