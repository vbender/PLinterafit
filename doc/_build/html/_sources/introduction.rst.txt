============
Introduction
============


"PLinterafit" stands for "PL spectra interactive fitting with Voigt functions". This Python package provides an interactive fitting tool for a Voigt function fitting to a photoluminscence spectrum 
where one can choose initial seed values of the Voigt function and analyse statistically the fitting quality of the optimized fit for the given seed values. 
The interactive slides help to refine the initial seed values which can lead to the best possible fitting results. 
Once the parameters of all Voigt functions are set, one can proceed to the next step of statistical analysis. 
Here, 'curve_fit' function from SCIPY package is used to obtain the best fit based on the seed values, and the corresponding goodness of fit values are estimated in terms of :math:`\chi^2` and p-value. 
At the end of the analysis, a figure will be produced including the measurement data, the fitted Voigt functions, and their total sum. The residue plot is also shown on the bottom of the figure
illustrating the difference between the measured value and the total fit.

This first version of the code was developed and used for the data analysis of a manuscript titled "**Spectroscopic study of the excitonic structure in monolayer** MoS :math:`_{2}` **under multi-variate 
physical and chemical stimuli**" authored by *Viktor Bender*, *Tobias Bucher*, *Mohammad Nasimuzzaman Mishuk*, *Yuxuan Xie*, *Isabelle Staude*, *Falk Eilenberger*, *Kurt Busch*, *Thomas Pertsch*, and *Bayarjargal N.Tugchin*
as a joint collaboration between Humboldt University of Berlin and Friedrich Schiller University Jena. 
The project was funded by the Deutsche Forschungsgemeinschaft (DFG, German Research Foundataion) - CRC/SFB 951 HIOS "Hybrid Inorganic/Organic Systems for Opto-Electronics" B10 - Project Number 182087777, 
the DFG - CRC/SFB 1375 NOA "Nonlinear Optics down to Atomic scales" (Projects A3, B2, B3, and C2), the DFG - Germany’s Excellence Strategy - EXC 2051 - Project Number 390713860.

