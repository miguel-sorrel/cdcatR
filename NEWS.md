cdcat 1.0.4:
* Some examples removed from gen.data()
* Fixed attribute labels in att.plot()
* Fixed a bug for LR.2step() where the function stopped if the computations resulted in NA

cdcat 1.0.3:
* cdcat.summary: Minor modifications in the plots
* gen.itembank: Two additional arguments incorporated (gs.parm and catprob.parm) to allow for more targeted manipulation of the quality of the item bank
* cdcat: New arguments included for specification of the starting rule (startK and startRule), item exposure control (itemExposurecontrol, b, and maxr), and the inclusion of additional constraints (itemConstraint and constraint.args$ATTRIBUTEc)
* helper: Functions updated to accommodate the above-mentioned modifications

cdcat 1.0.2:
* Fixed a bug for cdcat() when ItemSelect = "JSD"
* Fixed a bug for cdcat() when ItemSelect = "PWKL"
* Minor modification of cdcat() for gathering the information (Q-matrix)

cdcat 1.0.1:
* Addressed the noLD unit test issues in the previous version

cdcat 1.0.0:

DESCRIPTION:
* version and date updated
* P. Nájera listed as a co-author
* package dependencies updated

att.plot:
* CDM package accepted
* nonparametric CD-CAT accepted
* k argument included for selecting the attributes to be plotted. If k = NULL (by default), all attributes are plotted
* color has been added to the plot: red = non-mastered, green = mastered, blue = unclassified

cdcat: 
* CDM package accepted
* several item selection rules included: GDI, JSD, MPWKL, PWKL, NPS, random (argument itemSelect)
* nonparametric CD-CAT included (itemSelect = "NPS")
	* different arguments available: gate, pseudo.prob, w.type, and seed (argument NPS.args). pseudo.prob and w.type are experimental
	* some arguments needed to be modified / included (e.g., Q is required with nonparametric CD-CAT)
* parallelization included (argument n.cores)
* some warning / error messages included
* information regarding multiple modes (for ML and MAP) included in the output
* i.print is a deprecated argument. When parallelizing, no information is printed
* updated information and references
* new example for the nonparametric CD-CAT included
* examples updated 

cdcat.summary:
* major update
* item exposure information included
* optional label argument included (in case two objects share model)
* CDM package accepted
* nonparametric CD-CAT accepted
* plots improved (e.g., without secondary grids)
* violin plot for variable length CD-CAT improved (e.g., geom_dotplot instead of geom_jitter, x axis labels removed)

cdcat.comp:
* function removed. cdcat.comp is now included in cdcat.summary

gen.data:
* new function to simulate data

gen.itembank:
* major update
* data generation under different CDMs using the GDINA::simGDINA function
* additional specifications when creating a Q-matrix included 

LR2step:
* example updated
* CDM package accepted
* error when K^*j = 4 fixed
* p.adjust.methods method and alpha.level arguments included
* new output included LR2.adjp and models.adjp

sim180combination, 180DINA, sim180GDINA
* updated

helper:
* all functions required for cdcat and gen.itembank functions included

zzz:
* version updated
