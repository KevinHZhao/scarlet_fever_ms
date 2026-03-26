# Finding the PD bifurcation

1.  Unfortunately, these steps (PD and $\alpha$ continuation) will have to be repeated for each segment you are interested in analyzing.
2.  Copy paste the values from `avg_beta_parms.txt` into the relevant section of `SIRmap_PD.c` (according to the existing template), comment out all segments besides the one you wish to analyze, and then run `UbunCreateSO.sh` to generate the library for XPPAUT.
3.  Run `AutoSIR_PD.ode` with XPPAUT.
4.  We want to converge to a period 1 attractor, so $R_0$ should be relatively low (I used $R_0 = 7$).
5.  Go (I -\> G), then repeat from last (I -\> L, several times), until converged to period 1 attractor (Data window shows same rows for I and S).
6.  Auto (F -\> A), Run (R), look for the period doubling bifurcation (where red becomes black), can ABORT once we pass it.
7.  Grab (G), press TAB until we reach the PD bifurcation (Type should be HB), then ENTER to confirm we want to grab that point. Can also click near the PD bifurcation to lock on to it instead of TAB.
8.  Close AUTO, go back into XPP, File (F), Save pars (S), at the bottom of our saved .par file is the $R_0$ value where the PD bifurcation occurs. We save in XPP instead of AUTO since the file saved in AUTO seems to truncate several decimal places.
9.  This is the value of $R_0$ we want to use when trying to find $\alpha$ for the sinusoidal model.

# $\alpha$ continuation

8.  Edit `AutoSIR_alpha.ode` so that in the parameter values, R0 is set to the value of $R_0$ we found earlier. I would assume it shouldn't be necessary to use two distinct .ode files, but when I tried doing the steps below for $\alpha$ continuation after finding the PD bifurcation in the same AUTO window, they didn't work.
9.  Run `AutoSIR_alpha.ode` with XPPAUT.
10. The value of $\alpha$ was purposely set to $0.09$ instead of $0.1$, which is the value used in finding the PD bifurcation. This is so that we can use XPPAUT to converge to a period 1 attractor (I -\> G, then I -\> L several times until the Data window has repeating rows), then run AUTO (F -\> A), and Run (R) until we pass $\alpha = 0.1$, where there should be a bifurcation that we can Grab (G, click, Enter).
11. With this $\alpha$ grabbed, go into two par plot mode (A -\> T), set Main Parm: `p`, Secnd Parm: `alpha`, Xmin: `0`, Ymin: `0`, Xmax: `1`, Ymax `1`, then Ok.
12. Run (R), and the $\alpha$ continuation should be plotted.
13. Write the pts (F -\> W). Unfortunately, I couldn't find any obvious ways of finding $\alpha$ when $p = 1$ in the plot from AUTO, so I simply interpolated from the written points (in R).

# Generate macpan/sin brute force data

1.  Go into `bfd/` and run the `create_bfd.R` and `create_bfd_sin.R` scripts to generate several .ode files.
2.  WARNING: This step is very slow and could take more than an hour, there are definitely faster methods of generating such brute force data, see for example [Simon Frost's julia code](https://github.com/epirecipes/sir-julia/blob/master/markdown/ode_bifurcation_bruteforce/ode_bifurcation_bruteforce.md): Run `run_bfd.sh` and `run_bfd_sin.sh` in that same folder to generate the macpan brute force data.

# Creating a bifurcation diagram

-   Can follow the steps from Krylova & Earn and Papst & Earn (using `AutoSIR_PD.ode`, copy pasting the relevant icsets), some improvements:

    -   ABORT/[Esc] doesn't stop due to AUTO thinking there's a connection to another branch (I think this is something like a mistaken PD bifurcation). Solve by changing the EPS\* settings in Numerics in AUTO, I used 1e-12. (write up more clearly what this means)
    -   Can use -O3 flag when compiling the SIRmap shared object (see UbunCreateSO.sh), not sure if it makes it faster but I feel like it did.
