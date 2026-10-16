(Oct 13th 2026)
https://github.com/larseeri/stoch-MCFNDP-gener


## Generator of instances of the linear, two-stage, stochastic, multi-commodity, capacitated, fixed charge network design problem with recourse

## README

This is the companion repository of Larsen et al. (2023) containing the GenStoch generator of instances of the linear, two-stage, stochastic, multi-commodity, capacitated, fixed charge network design problem (MCFNDP) with recourse. GenStoch generates finite sets of realizations, i.e., scenarios, while matching target moments and correlations, for specified subsets of parameters of a given instance of the deterministic MCFNDP. Details of this process are governed by user-specified configuration settings. The code implementing the algorithm presented in Hoyland et al. (2003) and built by Michal Kaut and Diego Mathieu is used to compute the scenarios. For detailed background, see Larsen et al. (2023).

### Content

The repository contains the following:

/USERGUIDE_STOCH.md: User guide detailing installation, operation, functionalities and usage of options of the generator.

/TUTOREXAMPLE_STOCH.md: Tutorial example about the use of the generator.

/program: Folder structure and files required to build and run the generator from the Linux command line.

/docker: Folder structure and files required to run the generator through a Docker container built from the supplied Docker image.

/library: Folder structure containing a generic library of stochastic MCFNDP instances generated with the GenStoch generator, based on deterministic MCFNDP instances generated with the GenDeterm generator (1).

/library/DESCLIBRARY_STOCH.txt: Detailed description of the library.

/external_sources: Folder structure containing the material pertaining to the randomization of deterministic instances that are not originating from the GenDeterm generator.

/AUTHORS.txt: Authors of the generator.

/LICENSE.txt: Statement of the license covering the generator.

(1) https://github.com/larseeri/determ-MCFNDP-gener


### Learning

To become proficient with the use of the GenStoch stochastic generator, we suggest to start with steps 1-2 and to iterate on steps 3-5:

1 Read Sections 1 and 3 of Larsen et al. (2023) attentively.
2 Skim the user guide USERGUIDE_STOCH.md entirely to build a cognitive structure (know where to seek information, start locating the information related to the task at hand). The user guide details installation, operation, functionalities and available options of the generator.
3 Study the tutorial example TUTOREXAMPLE_STOCH.md in detail.
4 Study USERGUIDE_STOCH.md in detail. 
5 Launch generation experiments.


### Bibliography

Høyland, K., Kaut, M., Wallace, S.W., A Heuristic for Moment-matching Scenario Generation. Computational Optimization and Applications, 24 (2-3), pp. 169–185, 2003

Larsen, E., Bisaillon, S., Cordeau, J.F. and Frejinger, E. Pseudo-random Instance Generators in C++ for Deterministic and Stochastic Multi-commodity Network Design Problems. ArXiv e-prints arXiv:2312.10717, 2023.

