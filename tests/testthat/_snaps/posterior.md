# continuous predictive runs

    Code
      stats
    Output
      # A tibble: 3 x 5
         dose reference_dose  mean `2.50%` `97.50%`
        <dbl>          <dbl> <dbl>   <dbl>    <dbl>
      1     1              0  1.99   0.535     3.91
      2     3              0  6.11   3.98      8.12
      3     5              0 10.7    8.72     12.6 

# binary predictive runs

    Code
      stats
    Output
      # A tibble: 3 x 5
         dose reference_dose  mean `2.50%` `97.50%`
        <dbl>          <dbl> <dbl>   <dbl>    <dbl>
      1     1              0  0.18  0.0225    0.4  
      2     3              0  0.31  0.0450    0.578
      3     5              0  0.23  0.0225    0.478

