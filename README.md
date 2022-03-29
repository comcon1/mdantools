# MDAnTools

Crafty scripts for MDAnalysis expanding standard GROMACS tools.

1. **fistat_mda.py** - analogous of the 'gmx potential' classical utility. Finally, the electron density is integrated twice to get the potential.
2. **gdens_mda.py** - analogous of the 'gmx density' classical utility. 
3. **h2ord_mda.py** - analogous of the 'gmx h2order' but calculating both 1st and 2nd order parameters (Aman, 2003, BJ). The anlge is calculated betweeen OZ axis and the dipole vector of every water molecule. Note that here you must specify names of oxygen and hydrogens of water molecule in your system.
4. **./lipidnet** - intermediate steps for generation of triangulated network demonstrating lateral positioning of lipid molecules in a monolayer (Ermakov et al., 2019: Supplementary Video).

**NOTE:** The difference in **\_mda**-utils from corresponding GROMACS utils is the reference plane: here we attach the plane to the mean position of some lipid atom. The rest is the same: space is sliced with selected distance gap and a measure is calculated in every slice.

# Citation

1. Molotkovsky, R. J., Galimzyanov, T. R., Khomich, D. A., Nesterenko, A. M., & Ermakov, Y. A. (2021). Inhomogeneity of polylysine adsorption layers on lipid membranes revealed by theoretical analysis of electrokinetic data and molecular dynamics simulations. Bioelectrochemistry, 141, 107828. https://doi.org/https://doi.org/10.1016/j.bioelechem.2021.107828
2. Ermakov, Y. A., Asadchikov, V. E., Roschin, B. S., Volkov, Y. O., Khomich, D. A., Nesterenko, A. M., & Tikhonov, A. M. (2019). Comprehensive Study of the Liquid Expanded–Liquid Condensed Phase Transition in 1,2-Dimyristoyl- sn -glycero-3-phospho- l -serine Monolayers: Surface Pressure, Volta Potential, X-ray Reflectivity, and Molecular Dynamics Modeling. Langmuir, 35(38), 12326–12338. https://doi.org/10.1021/acs.langmuir.9b01450
