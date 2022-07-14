# print methods

    Code
      print(mod)
    Output
      -------------
      dreamer Model
      -------------
      type    : linear
      response: continuous
      
      Dose Response
      |hyperparameter | value|
      |:--------------|-----:|
      |mu_b1          | 0.000|
      |sigma_b1       | 1.000|
      |mu_b2          | 0.000|
      |sigma_b2       | 1.000|
      |shape          | 1.000|
      |rate           | 0.001|
      |w_prior        | 1.000|

---

    Code
      print(output)
    Output
      -----------------
      dreamer Model Fit
      -----------------
      doses: 0, 2.5, 5
      
      | dose|   mean|  2.50%| 97.50%|
      |----:|------:|------:|------:|
      |  0.0|  0.177| -0.296|  0.819|
      |  2.5|  5.593|  5.281|  5.894|
      |  5.0| 11.009| 10.529| 11.405|
      
      |model      | prior weight| posterior weight|
      |:----------|------------:|----------------:|
      |mod_linear |            1|                1|

---

    Code
      print(output$mod_linear)
    Output
      -----------------
      dreamer Model Fit
      -----------------
      doses: 0, 2.5, 5
      
      | dose|   mean|  2.50%| 97.50%|
      |----:|------:|------:|------:|
      |  0.0|  0.177| -0.296|  0.819|
      |  2.5|  5.593|  5.281|  5.894|
      |  5.0| 11.009| 10.529| 11.405|

---

    Code
      print(mod_binary)
    Output
      -------------
      dreamer Model
      -------------
      type    : linear
      response: binary
      
      Dose Response
      |hyperparameter |value  |
      |:--------------|:------|
      |mu_b1          |0      |
      |sigma_b1       |1      |
      |mu_b2          |0      |
      |sigma_b2       |1      |
      |w_prior        |1      |
      |link           |probit |

