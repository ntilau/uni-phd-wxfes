rem fe WG_p1 1e10 ++ +dbl -verbose +hfss +msh
rem fe WG_p1 1e10 ++ +dbl -verbose
rem fe WG_p1 1e10 ++ +tfe -verbose
rem fe WG_p1 1e10 +dd 8 +gmres 1e-4 100 ++ +dbl -verbose
rem fe WG_p1 1e10 +dds 8 +tfe +gmres 1e-4 100 ++ +dbl -verbose

rem fe WG_p1 1e10 ++ +dbl -verbose +hfss +href qa0.0000000001 +msh
rem fe WG_p1 1e10 ++ +dbl -verbose
rem fe WG_p1 1e10 ++ +tfe -verbose
fe WG_p1 1e10 +dd 14 +gmres 1e-4 100 ++ +dbl -verbose
fe WG_p1 1e10 +dd 14 +tfe +gmres 1e-4 100 ++ +dbl -verbose

fe WG_p1 1e10 ++ +dbl -verbose +hfss +href qa0.00000000001 +msh
fe WG_p1 1e10 ++ +dbl -verbose
fe WG_p1 1e10 ++ +tfe -verbose
fe WG_p1 1e10 +dd 18 +gmres 1e-4 100 ++ +dbl -verbose
fe WG_p1 1e10 +dd 18 +tfe +gmres 1e-4 100 ++ +dbl -verbose
