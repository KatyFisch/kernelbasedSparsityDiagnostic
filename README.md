# kernelbasedSparsityDiagnostic
This repository contains the code to the simulation studies of the paper "A Diagnostic to Find and Help Combat Positivity Issues – with a Focus on Continuous Treatments" by Katharina Ring and Michael Schomaker

The analysis is conducted in ```analysis.R```, which simulates data and derives the true curves and estimated curves for each estimand using ```curvecalc_functions.R```, plots them using ```plot_functions.R``` and assembles the plots into pdfs in ```pdf_functions_wdiagnostic.R```. For the diagnostic, the file ```diagnostic_function.R``` is used.
