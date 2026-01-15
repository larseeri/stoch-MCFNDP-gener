// Copyright (c) 2025 Eric Larsen
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the “Software”), to
// deal in the Software without restriction, including without limitation the
// rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
// sell copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
// IN THE SOFTWARE.
// ****************************************************************************
// ****************************************************************************
//
// SCENARIOS GENERATOR FOR THE LINEAR, TWO-STAGE, MULTI-COMMODITY, FIXED
// CHARGE, CAPACITATED NETWORK DESIGN PROBLEM WITH STOCHASTIC RECOURSE
//
// This application randomly generates finite sets of random realizations,
// i.e., scenarios, while matching target moments and correlations, for
// specified subsets of parameters of a given instance of the deterministic
// multi-commodity, fixed charge, capacitated network design (MCFND) problem.
// Details of this process are governed by user-specified configuration
// settings. The code implementing the algorithm presented in Hoyland, Kaut
// and Wallace (2003) (HKW) and built by Michal Kaut and Diego Mathieu is used
// to compute the scenarios.
//
// For detailed background, see
// Larsen, E., Bisaillon, S., Cordeau, J.F. and Frejinger, E.
// Pseudo-random Instance Generators in C++ for Deterministic and Stochastic
// Multi-commodity Network Design Problems.
//
// The two following paragraphs are partly adapted from
// ../HOYLAND_KAUT_WALLACE/sg_HKW.c:
//
// Scenario generation code is based on 'A Heuristic for Moment-matching
// Scenario Generation' by K. Høyland, M. Kaut & S.W. Wallace; Computational
// Optimization and Applications, 24 (2-3), pp. 169–185, 2003;
// doi:10.1023/A:1021853807313. Code by Michal Kaut (michal.kaut@iot.ntnu.no)
// & Diego Mathieu.
//
// The HKW code generates scenarios for multivariate random variables by
// matching specified first four moments and correlations. The algorithm
// is iterative and gradually closes the gap from target values. Although
// convergence has not been demonstrated theoretically, extensive empirical
// testing and application have shown that convergence occurs provided the
// requested number of scenarios is sufficiently large. Several trials may be
// required, each one starting from new random initializations automatically.
// (Cf. Section 6 below.)
//
// 1 Building and running the program directly from the Linux command line
// ***********************************************************************
//
// 1.1 Building the program from the Linux command line
// ****************************************************
//
// Under Linux, follow these steps to build the program:
//
// a- Install the CPLEX solver and/or the HIGHS solver (if not already done)
// and note the installation folder(s). According to our verifications, CPLEX
// and HIGHS work equally well in our application. HIGHS is freely available
// under the MIT license at https://highs.dev/. See Huangfu, Q., Hall, J.A.J.
// Parallelizing the dual revised simplex method. Math. Prog. Comp. 10, 119-142
// (2018).
// b- Either (i) use Git to clone locally the remote repository
// https://github.com/larseeri/stoch-MCFNDP-gener, or (ii) generate, download
// locally to an appropriate folder and uncompress a .zip archive containing
// the remote repository file tree.
// c- Edit the file named "makefile" to set the variables indicating the
// locations relevant to the CPLEX solver and/or the HIGHS solver, as noted in
// step a-. These locations are chiefly CPLEX_HOME and HIGHS_HOME. Save the
// modifications.
// d- In a Linux terminal, go to the /program folder just created.
// e- Make sure that the gcc C++ compiler and the make C++ build tool are
// available.
// f- In the terminal, run ONLY ONE of these commands to build the program
// equipped either with the CPLEX solver or with the HIGHS solver:
//
// make CPLEX_SOLVER=1
//
// make HIGHS_SOLVER=1
//
// An executable file named "stoch_gen" will be created inside the /program
// folder. Building with versions 11.5.0, 12.2.0, 13.2.0, 14.2.0 and 15.1.0 of
// the gcc compiler has been verified.
//
// Remarks:
// - To build once again with HIGHS if the program has already been built with
// CPLEX, run the following:
//
// make clean CPLEX_SOLVER=1
// make HIGHS_SOLVER=1
//
// - To build once again with CPLEX if the program has already been built with
// HIGHS, run the following:
//
// make clean HIGHS_SOLVER=1
// make CPLEX_SOLVER=1
//
// 1.2 Running the program with Linux command line instructions
// ************************************************************
//
// The program can be run as follows from the /program folder just created
// along with default option parameter values:
//
// ./stoch_gen
//
// Options for running the program can be specified directly in the command
// line or in a configuration file. By default, the program will attempt to
// read option values from file stoch_config.par. (This default setting can be
// updated by modifying line 1197 of file main.cpp.) The default file can
// be overridden with option +F:
//
// ./stoch_gen +F newParamfile.par
//
// This assumes that stoch_config.par and newParamfile.par are located at the
// root of the /program folder. Otherwise, their absolute paths or their paths
// relative to the /program folder must be explicited.
//
// The default configuration file stoch_config.par has initially been filled
// with illustrative values and should be edited appropriately if no
// replacement configuration file is indicated on the command line with option
// key +F.
//
// Verbosity option +v should be used to obtain more information about all
// options that have been set in the parameter file or in the command line:
//
// ./stoch_gen +v +F newParamfile.par
//
// The configuration file should have this format:
//
// optionKey1   value1
// optionKey2   value2
// .
// .
// optionKeyN   valueN
//
// It may specify at most one option per line. Option keys may be prefixed or
// not with '-'. Comments may appear either at the beginning of a line, in
// which case they must start with a non-alphabetic character (but not with
// '-'), or at the end of a line.
//
// The general command line for running the program is as follows:
//
// ./stoch_gen [+v] [+F nameOfConfigFile] [main and HKW options]
//
// Once the configuration file has been processed, options appearing in the
// command line (other than +F and +v) are processed. Option values specified
// in the command line supersede those specified in the configuration file.
// That is, if included, they replace, modify and supplement the latter.
// Option keys appearing in the command line must be prefixed with '-' as
// follow:
//
// -optionkeyA valueA -optionKeyB valueB...
//
// By default, boolean values are set to false and presence of option key is
// sufficient to set the corresponding value to true. For other options
// expecting string, integer or double values, values should be specified.
//
// To display the complete documentation about the application and the
// available options, run
//
// ./stoch_gen -H
// or
// ./stoch_gen -h
//
// Summing up, option values are successively updated from these locations in
// this order:
// (i) hard-coded initialization values in lines 27-51 of configuration.cpp
// (ii) values appearing in configuration file whose name is specified in
// command line with +F (if it is specified and if it exists); if no valid
// configuration file is specified with +F in command line, values will be
// taken from default configuration file whose name is specified in line 1197
// of main.cpp; at this point, the latter should exist
// (iii) values specified in the command line.
// Notice that location (i) is insufficiently transparent to the user and
// should not be relied upon. Parameter values should be set through (ii) and 
// (iii).
//
// Remark:
// See Section 8 below for examples of configuration files and instructions
// when the program is run with command line instructions.
//
// 2 Running the program through Docker instructions
// *************************************************
//
// An advantage of running the executable stoch_gen inside a Docker container
// is the possibility to avoid selection and installation of suitable C++
// compiler, C++ make build tool and required C++ libraries. It suffices that
// the Docker engine be installed on the host where the computation will be
// run. For this, see https://docs.docker.com/engine/install/.
//
// Once the Docker engine is installed, follow these steps:
//
// a- Open a terminal and retrieve from the GitHub container registry the
// Docker image ghcr.io/larseeri/stoch_gen_image:latest as follows:
//
// docker pull ghcr.io/larseeri/stoch_gen_image:latest
//
// (From this image, Docker containers running the executable 'stoch_gen' will
// be generated.)
// b- For simplicity, define the alias 'stoch_gen_image' through
//
// docker tag ghcr.io/larseeri/stoch_gen_image:latest stoch_gen_image
//
// Remark:
// If the Docker image stoch_gen_image were instead available as the tar ball
// stoch_gen_image.tar, then it could be loaded through
//
// docker load --input stoch_gen_image.tar
//
// c- Either (i) use Git to clone locally the remote repository
// https://github.com/larseeri/stoch-MCFNDP-gener, or (ii) generate, download
// locally to an appropriate folder and uncompress a .zip archive containing
// the remote repository file tree. The /docker folder contains (i) the
// primitive files named 'Dockerfile' and 'makefile' based on which the
// Docker image was originally generated, (ii) a subfolder /inout used
// for exchanging files between the program running inside a Docker container
// and the directory structure of the host, (iii) some auxiliary files, and
// (iv) the subfolder /examples containing a number of configuration files and
// their resulting outputs (examined in Section 5 below).
// d- In a terminal, go to the /docker folder just created.
//
// Then, instructions similar to those appearing in Section 1.2 can be run from
// the terminal on AMD and Intel CPUs (other CPUs are addressed further below).
// For example, the default Docker instruction
//
// docker run --rm --mount type=bind,src=absolutePathToDockerFolder/inout,dst=/inout stoch_gen_image /stoch_gen
//
// is equivalent to the default Linux command line instruction
//
// ./stoch_gen
//
// In general, Docker instructions can be formulated from those appearing in
// Section 1.2 by (i) adding the prefix 
//
// docker run --rm --mount type=bind,src=absolutePathToDockerFolder/inout,dst=/inout stoch_gen_image
//
// and by (ii) replacing the executable location ./stoch_gen with /stoch_gen
// (i.e. removing the initial dot), since the executable stoch_gen is located
// at the root of the Docker container being generated and since locations 
// inside a Docker container must be accessed through absolute paths.
//
// Hence, the following command line instructions:
//
// ./stoch_gen +v +F newParamfile.par
//
// ./stoch_gen -h
//
// ./stoch_gen +v +F newParamFile.par -r 1234 -R 4567
//
// respectively correspond to the following Docker instructions:
//
// docker run --rm --mount type=bind,src=absolutePathToDockerFolder/inout,dst=/inout stoch_gen_image \
// /stoch_gen +v +F /inout/newParamfile.par
//
// docker run --rm --mount type=bind,src=absolutePathToDockerFolder/inout,dst=/inout stoch_gen_image \
// /stoch_gen -h
//
// docker run --rm --mount type=bind,src=absolutePathToDockerFolder/inout,dst=/inout stoch_gen_image \
// /stoch_gen +v +F /inout/newParamFile.par -r 1234 -R 4567
//
// Important: 
// When running the program through Docker instructions, all input and output
// files must be located inside the docker/inout subfolder. Within the Docker 
// instructions and the configuration files, this subfolder must be denoted
// /inout since it is located at the root of the Docker container. Hence, a
// replacement configuration file located at docker/inout/newParamFile.par
// will be specified in the Docker instruction with +F /inout/newParamFile.par.
// Similarly, a base network file at docker/inout/baseNetworkFile.std will be
// specified in the Docker instruction or inside the configuration file with
// -I /inout/baseNetworkFile.std.
//
// The default configuration file is located at docker/inout/stoch_config.par.
// It has initially been filled with illustrative values and should be edited
// appropriately if no replacement configuration file is indicated on the
// Docker instruction line with option key +F.
// 
// Option parameter values are successively updated from these locations in
// this order:
// (i) initialization values hard-coded in lines 27-51 of configuration.cpp,
// (ii) either replacement configuration file specified in the Docker
// instruction line with option key +F, or default configuration file 
// ./inout/stoch_config.par if no replacement configuration file is 
// specified in the Docker instruction line with option key +F,
// (iii) values specified in the Docker instruction.
// Notice that location (i) is insufficiently transparent to the user and
// should not be relied upon. Parameter values should be set through (ii) and
// (iii).
//
// Important:
// The supplied Docker image conforms to the amd64 CPU instructions standard
// and will run natively on AMD and Intel CPUs. On CPUs conforming instead to
// the arm64 standard, e.g., on recent Apple computers, the supplied Docker 
// image can be run through emulation. Essentially, this may be achieved as
// follows:
// a- Running the following preliminary instruction once. (For details, see
// https://hub.docker.com/r/tonistiigi/binfmt.)
//
// docker run --privileged --rm tonistiigi/binfmt --install all 
//
// b- Adding the option --platform linux/amd64 to the Docker instruction line
// prefix as follows:
//
// docker run --platform linux/amd64 --rm --mount type=bind,src=absolutePathToDockerFolder/inout,dst=/inout stoch_gen_image
//
// The following characteristics are specific to the operation of the program
// through Docker instructions:
// - An empty folder absolutePathToDockerFolder/inout must already exist on
// the host for the Docker instructions to work correctly.
// - Option --mount type=bind,src=absolutePathToDockerFolder/inout,dst=/inout
// connects subfolder docker/inout on the host to folder /inout
// of the Docker container. Hence, the connection between
// absolutePathToDockerFolder/inout on the host and folder /inout in the
// container is used to supply input files to and retrieve output files from
// the container.
// - Option --rm ensures suppression of the container implicitly generated by
// docker run. Note that the container created by docker run for running the
// executable stoch_gen stops running immediately after the latter finishes
// and that it cannot be reused afterwards, for instance through a docker exec
// instruction. Immediate suppression with --rm avoids uselessly cluttering
// Docker's container repository.
//
// The following Docker instructions may be useful:
//
// docker image ls (list all currently existing docker images)
//
// docker image rm --force imageName (forcefully remove imageName)
//
// docker container ls (list all running containers)
//
// docker container ls --all (list all containers, running or not)
//
// See Section 8 below for examples of configuration files and instructions
// when the program is run through Docker.
//
// 3 Randomized elements, target moments and target correlations
// *************************************************************
//
// The generator rests on three primary components: a subset of randomized
// parameters, target moments and target correlations. We proceed to describe
// how each one is specified.
//
// Remarks:
// - Section 7 below describes detailed usage of every option. Section 8
// presents systematically some examples of configuration options specified in
// configuration files.
// - For the purposes of the current discussion, <...> indicates that some
// value or values of the type specified inside the brackets is expected. For
// instance, <string> means that a string is expected whereas <> means that no
// value is expected.
//
// 3.1 Base deterministic MCFND instance and randomized elements
// *************************************************************
//
// Option S <int> must always be specified, either in the instruction line or
// in the configuration file. It indicates which subset of parameters of a
// given base deterministic MCFND instance will be assigned a joint probability
// distribution and will thus vary between scenarios. We call these parameters
// randomized elements. Option S expects an integer from 1 to 31 that is the
// sum of the following, where symbol -> is a shorthand notation for "calls
// for randomization of":
// 1 -> demand volumes
// 2 -> total capacities of arcs
// 4 -> commodity-specific capacities of arcs
// 8 -> fixed costs
// 16 -> variable costs (either commodity-specific or not)
// Hence, specifying -S 3 in the instruction line indicates that total
// capacities of arcs (2) and demand volumes (1) must be randomized
// (2 + 1 = 3).
//
// By default, the name of the file describing the base deterministic MCFND
// network is './base_network.dow' (Linux instructions),
// './inout/base_network.dow' (Docker instructions) and the default input
// format is DOW. This name may be changed with option I <string> and the input
// format can be indicated accordingly with option F <char>. For example,
// specifying -S 3 -I instB.std -F G in the instruction line indicates that
// demands and total capacities of arcs will be randomized based on the
// deterministic network stored in file instB.std using the generic STD format.
// (Cf. Section 5.1 below for details.)
//
// Important:
// It is critical that the content of the file describing the base
// deterministic MCFND instance and the indicated format of its content match
// together.
//
// Remark:
// While the generator is sufficiently general to randomize fixed costs, the
// latter are viewed as being non-stochastic in the context of two-stage
// stochastic programming. Hence, in standard applications, randomization of
// the MCFND will be limited to a subset of the following: demand volumes,
// total capacities of arcs, commodity-specific capacities of arcs and variable
// costs (either commodity-specific or not).
//
// 3.2 Target moments
// ******************
//
// Target moments can be set and saved as follows.
// a- If option key G is included, target moments will be generated according
// to a distributional assumption and will be written individually in a file
// whose default name is './target_moms.dat' (Linux instructions) or
// './inout/target_moms.dat' (Docker instructions). This name can be specified
// otherwise with option key M. The distributional assumption can be specified
// with option key T <char> where <char> stands for either U (uniform) or D
// (triangular) and two parameters, alpha and beta, that can be specified using
// options A <double> and B <double>.
//
// b- If option key G is not included, individual target moments will be read
// from the file whose default name is is './target_moms.dat' (Linux
// instructions) or './inout/target_moms.dat' (Docker instructions). The latter
// can be specified otherwise with option key M. Format of the target moments
// file is specified in Section 5.2.1 below. Without option G, option T will be
// rejected and options A and B will be ignored.
//
// Important:
// In case b-, that is when option G is omitted and target moments are read
// from a file, it is crucial that they are compatible with a proper
// probability distribution. If this condition fails, an attempt at
// randomization is meaningless and convergence of the HKW algorithm is
// expected to fail. In case a-, that is when option G is included, the target
// moments that are generated are by construction compatible with a proper
// probability distribution.
//
// 3.3 Target correlations
// ***********************
//
// Target correlations can be set and saved as follows.
// a- If option key K is included, target correlations can be directly
// specified in blocks either (i) in the instruction line, or (ii) in the
// configuration file. A common target correlation coefficient within subset of
// parameters U or between subsets of parameters U and V is specified with
// X<UU> <double> and X<UV> <double> respectively, where <double> is the target
// correlation coefficient and U and V, can take these values, where symbol ->
// is a shorthand notation for "stands for":
// D -> demands
// A -> total capacities of arcs
// C -> commodity-specific capacities of arcs
// F -> fixed costs
// V -> variable costs
// Target correlations that are left unspecified are assumed by default to be
// zero. The resulting individual correlations will be written to a file whose
// default name is './target_corrs.dat' (Linux instructions) or
// './inout/target_corrs.dat' (Docker instructions). This name can be
// specified otherwise with option key C.
//
// Important:
// Cases (i) and (ii) are not coexistent: if any correlation information is
// supplied in the instruction line, any correlation information supplied in
// the configuration file is ignored.
//
// b- If option key K is not included, individual correlations will be read
// from a file whose default name is './target_corrs.dat' (Linux instructions)
// or './inout/target_corrs.dat' (Docker instructions). The latter can be
// specified otherwise with option key C. Format of the target correlations
// file is specified in Section 5.2.2 below. Without option K, options X?? will
// be rejected.
//
// Important:
// In case a-, when option key K is included, it is crucial that the block
// correlations that are specified result in a positive definite correlation
// matrix. Similarly, in case b-, when option key K is omitted, the supplied
// matrix must define a positive definite correlation matrix. Whenever the
// alleged correlation matrix, either constructed (case a-) or read from file
// (case b-), fails positive definitess, its Cholesky decomposition will fail
// in the HKW algorithm. A suggestion to verify and correct the alleged
// correlation matrix will be displayed and the program will be interrupted.
//
// 3.4 Example
// ***********
//
// Consider this hypothetical Linux command line instruction:
//
// ./stoch_gen -S 3 -I instB.std -F G -K -XDD 0.5 -XDA -0.3 -XAA 0.7 -G -T U -A 0.25 -B 0.25
//
// It is interpreted as follows:
// -S 3: demands and total capacities of arcs will be randomized
// -I instB.std: base deterministic network is to be read from file instB.std
// -F G: format of the latter is generic STD
// -K: target correlations are specified in blocks (based on explicit and
// default values)
// -XDD 0.5: target correlations among demands are 0.5 for all dissimilar
// pairs and self-correlations are 1.0
// -XDA -0.3: target correlations between demands and total capacities of
// arcs are -0.3 for all pairs
// -XAA 0.7: target correlations among total capacities of arcs are 0.7
// for all dissimilar pairs and self-correlations are 1.0
// -G: target moments are generated (based on explicit and default
// specifications)
// -T U: target moments are generated under the assumption that
// randomizations each obey a uniform distribution
// -A 0.25: parameter alpha defining the latter distribution is 0.25
// -B 0.25: parameter beta defining the latter distribution is 0.25
//
// 4 Output file
// *************
//
// By default, the name of the file containing the MCFND instances generated
// by the application is './generated_networks.txt' (Linux instructions) or
// './inout/generated_networks.txt' (Docker instructions). This can be
// changed with option O. This text file superposes representations of complete
// deterministic MCFND instances, one for each scenario realization. Each
// representation is preceded by a separator as follows:
//
// SCENARIO   <int>
//
// Immediately after the separator, the generated MCFND instance corresponding
// to the current scenario is described in its entirety with a format (DOW,
// restricted STD or generic STD; cf. Section 5.1 below), matching that used to
// supply the base deterministic MCFND instance. Stacking and separating full
// descriptions of the MCFND realizations in this fashion makes it easy to
// pinpoint and read the first and second stage information relevant to every
// individual scenario.
//
// Remark:
// Feasibility of the second stage of each generated MCFND instance is verified
// with the selected solver (CPLEX or HIGHS) after opening all arcs. Only the
// generated networks that are demonstrated feasible are added to the output
// file. The application produces a report indicating how many scenarios were
// generated and how many among the latter are feasible.
//
// 5 Formats of data files
// ***********************
//
// 5.1 Base deterministic MCFND instance and generated MCFND instances
// *******************************************************************
//
// DOW and STD formats are commonly used to store deterministic MCFND instances
// in text files. These formats are described in files DOW-format-desc.txt and
// STD-format-desc.txt. The application can read from these two formats and, in
// addition, from a restricted version of STD. Thus, three formats are
// available for supplying the base deterministic MCFND instance: DOW, STD
// (calling it generic STD, or GSTD) and a restricted version of STD (calling
// it restricted STD, or RSTD). One of these three formats may be selected with
// option F <char>.
// a- DOW presupposes that (i) for each arc, capacities and variable
// costs are identical for all commodities and (ii) each commodity has a
// single source node (at which volume > 0) and a single sink node (at which
// volume < 0).
// b- Generic STD imposes neither (i), nor (ii).
// c- Restricted STD imposes only (ii).
//
// Important:
// When the base deterministic MCFND instance is received in the generic STD
// format, there are potentially multiple source and/or sink nodes for some or
// all commodities. In this case, whenever randomization of commodity volumes
// is requested through option S, it is handled through these steps:
// a- For each commodity, the total volume, the flows out of each source node
// and the flows into each sink node occurring in the base deterministic MCFND
// network are saved for ulterior computations while reading (see
// InstanceGStd::read() in generate.cpp).
// b- Total commodity volumes are randomized (along with other stochastic
// elements if there are any) as they are when there are single source and sink
// nodes for each commodity.
// c- For each commodity, the randomized flows out of each source node and
// randomized flows into each sink node are calculated by distributing the
// randomized total commodity volume using the original ratios of inflows and
// outflows over total volume obtained in step a-.
//
// Remark:
// The MCFND instances generated by the application (one for each scenario) are
// expressed in the same format as that used for supplying the base
// deterministic MCFND instance. Their descriptions are stacked and returned in
// the output file whose default name './generated_networks.txt' (Linux
// instructions), './inout/generated_networks.txt' (Docker instructions) can be
// modified with option O. (Cf. Section 4 above.)
//
// 5.2 Input files supplied to HKW algorithm
// *****************************************
// (partly adapted from ../HOYLAND_KAUT_WALLACE/sg_HKW.c)
//
// 5.2.1 Target moments
// ********************
//
// When option G is omitted from instruction line and from configuration file,
// target moments are read from a file named './target_moms.dat' (Linux
// instructions) or './inout/target_moms.dat' (Docker instructions) unless this
// is overridden with option M. The file containing the target moments matrix
// is formatted as follows (see Section 8 for examples):
// a- 1st line shows integer 4 (as there are 4 target moments per randomized
// stochastic element and therefore 4 rows in the matrix of target moments),
// b- 2nd line states the total number of individual randomized stochastic
// elements,
// c- 3rd line states the first target moment of each randomized stochastic
// element, separated by spaces and listed in this order: demands (if
// randomized), total capacities of arcs (if randomized), commodity-specific
// capacities of arcs(if randomized), fixed costs (if randomized), variable
// costs (if randomized). The number of figures in this line is equal to the
// number specified in 2nd line.
// d- 4th to 6th lines respectively state second, third and fourth target
// moments, exactly as 3rd row does for first target moment.
// e- Remaining lines of the file are ignored.
//
// Important:
// Exact definitions of the four target moments must satisfy these
// conditions:
// a- Population estimators are used (as in spreadsheets).
// b- 2nd moment is defined as variance, that is second central
// statistical moment, and not standard deviation.
// c- 3rd moment is defined as third central statistical moment divided
// by third power of standard deviation (i.e., standardized skewness).
// d- 4th moment is defined as fourth central statistical moment divided
// by fourth power of standard deviation minus 3 (i.e., excess standardized
// kurtosis).
//
// 5.2.2 Target correlations
// *************************
//
// When option K is omitted from instruction line and from configuration file,
// target correlations are read from a file named './target_corrs.dat' (Linux
// instructions) or './inout/target_corrs.dat' (Docker instructions), unless
// this is overridden with option C. The file containing the target
// correlations matrix is formatted as follows (see Section 8 for examples):
// a- 1st and 2nd lines are identical and show the total number of individual
// randomized elements.
// b- Remaining lines of the file contain the target correlations matrix, one
// row of the matrix per line of the file. Individual correlations are
// separated by spaces.
// c- The correlation pairs appearing in the matrix are ordered by row and by
// column. Order in a row and order in a column is similar to that in a row
// of the matrix of target moments (cf. Section 4.1.1).
// d- Any additional line in the file is ignored.
//
// 5.2.3 Probabilities of scenarios
// ********************************
//
// If option p is omitted, then scenarios are assumed to be equiprobable.
// Option p can be used to supply the name of a file containing probabilities
// of all scenarios. The following format must be satisfied (see supplied
// example probs_100_eq_example.dat):
// a- 1st line states the number of scenarios
// b- 2nd line states all probabilities, separated by spaces
// c- Remaining lines of the file are ignored.
// This file would be supplied as './probs_100_eq_example.dat' (Linux
// instructions) or './inout/probs_100_eq_example.dat' (Docker instructions).
//
// 5.2.4 Raw scenarios in HKW format
// *********************************
//
// Raw scenarios are saved in HKW format to './raw_scens.dat' (Linux
// instructions) or './inout/raw_scens.dat' (Docker instructions), unless this
// is overriden with option o. Raw scenarios saved in previous run can be used
// as starting values in a following run (using option s). These will only be
// used in the first trial of the current run of the HKW algorithm. If this
// trial does not converge and the maximum number of trials specified by option
// t is larger than 1, the next trial starts from random guesses.
//
// Important:
// Raw scenarios pertain to the individual stochastic elements. They are
// distinct from the MCFND instance scenarios produced by the application and
// described in Section 4 above.
//
// 6 Convergence
// *************
//
// Remark:
// The level of on-screen reporting by the HKW algorithm should be set at
// no less than the default value of 3 (option l) to adequately follow the
// convergence process.
//
// The HKW algorithm might fail to satisfy the specified target moments and
// correlations within the default tolerances with the requested number of
// scenarios within the default number of trials and iterations. If
// convergence fails, the following course of action is recommended:
//
// 1- Verify that target moments and correlations are proper.
//
// If option G has not been selected, target moments are read from a file. It
// should be ascertained they are compatible with a probability distribution.
// If this condition fails, convergence of the HKW algorithm is expected to
// fail. By construction, target moments that are generated when option G is
// selected are correctly specified.
//
// If option K has not been selected, a target correlation matrix is read from
// a file. it should be ascertained this matrix is well-specified (square,
// symmetric, all elements in [-1,1], diagonal elements equal to 1, positive
// definite). If this condition fails, the HKW algorithm is expected to fail.
//
// If option K has been selected, a target correlation matrix is built from
// block correlations specified with options X??. it should be ascertained the
// resulting alleged target correlation matrix is well-specified.
//
// Remarks:
// - If the alleged correlation matrix fails positive definiteness, the HKW
// algorithm initially (i.e. before any iteration) reports a failure of its
// Cholesky decomposition and displays an explanation.
// - However, failure of the Cholesky decomposition during the iterations of
// the HKW algorithm might be due to an insufficient requested number of
// scenarios (see Paragraph 3 below).
//
// 2- Verify that slow convergence is not responsible.
//
// If both target moments and target correlations are known to be properly
// specified, consider launching a new run after increasing the maximum numbers
// of trials and iterations, respectively options t and i. Doubling their
// values from default is suggested.
//
// Remark:
// Default values for the numbers of trials and iterations are usually
// sufficient to reach convergence under the default error tolerances. Failure
// of convergence is most frequently due to the causes examined in Paragraphs 1
// and 3 above and below.
//
// 3- Adjust the number of scenarios.
//
// If both target moments and target correlations are properly specified and
// slow convergence is excluded, then the requested number of scenarios is too
// small to support the specified target moments and target correlations.
// Hence, the requested number of scenarios, option n, may be progressively
// increased by a moderate factor (e.g., between 1.5 and 5) until the algorithm
// converges. Once it has converged, the number of requested scenarios may be
// progressively adjusted downward, if desired, while the algorithm is still
// converging.
//
// Remark:
// The lower bound on the number of requested scenarios leading to convergence
// is positively related to the following: number of random elements, magnitude
// of 3rd and 4th target moments, heterogeneity of target correlations,
// heterogeneity of probabilities.
//
// 7 Usage of options
// ******************
//
// (Some parts are adapted from ../HOYLAND_KAUT_WALLACE/sg_HKW.c.)
//
// Except for option key S (and leaving aside option key F when it is prefixed
// by '+' in the instruction line), non-boolean user-selected options that are
// specified neither in the instruction line nor the configuration file are set
// to the values that are hard-coded in the constructor in lines 27-51 of
// configuration.cpp. Boolean user-selected options that are specified neither
// in the instruction line nor in the configuration file are set to false by
// default.
//
// Notation:
// For the purposes of this presentation, <...> indicates that some value or
// values of the type specified inside the brackets is expected. For instance,
// <string> means that a string is expected whereas <> means that no value is
// expected.
//
// 7.1 Options available in instruction line only (Linux or Docker)
// ****************************************************************
//
// +F <string>; name of configuration file supplying option parameter values;
// this replaces the default file whose name is hard-coded in line 1197 of
// main.cpp; any option value supplied in instruction line replaces value of
// same option supplied in configuration file
//
// +v <> (no value supplied); requests verbosity while parsing options in
// instruction line or configuration file
//
// Remark:
// Unknown option keys prefixed with '+' in the instruction line are
// disregarded.
//
// 7.2 Options available in instruction line (Linux or Docker) and in configuration file
// *************************************************************************************
//
// The following options must be prefixed by '-' when they are called from the
// instruction line. This prefix can be omitted when the options are called
// from a configuration file. To supply an option, a line of a configuration
// file must begin with a dash '-' or a letter. A line beginning with any other
// character will be considered a comment. The remainder of a line beyond the
// option key and its value (if any is expected) is disregarded and may
// therefore be used for adding comments.
//
// Remark:
// Unknown option keys following '-' in the instruction line, stray characters
// or strings appearing without the prefix '-' in the instruction line and
// unknown option keys (following '-' or not) at beginning of line in a
// configuration file are disregarded.
//
// 7.2.1 Main generator options
// ****************************
//
// H <> (no value supplied); display complete documentation about the
// application and the available options
//
// h <> (no value supplied); display complete documentation about the
// application and the available options
//
// F <char>; format of file where deterministic MCFND base network is to be
// read; D: DOW, G: generic STD, R: restricted STD; default: D
//
// I <string>; name of file where deterministic MCFND base network is to be
// read; default value is './base_network.dow' (Linux instructions) or
// './inout/base_network.dow' (Docker instructions)
//
// O <string>; name of file where MCFND instances resulting from
// scenario generation are to be written; default value is
// './generated_networks.txt' (Linux instructions) or
// './inout/generated_networks.txt' (Docker instructions)
//
// S <int>; identifies which subsets of parameters should vary between
// scenarios; expects a number (from 1 to 31) which is a sum of the following:
// 1 -> demands
// 2 -> total capacities of arcs
// 4 -> commodity-specific capacities of arcs
// 8 -> fixed costs
// 16 -> variable costs
//
// K <> (no value supplied); if included, target correlations are directly
// specified in blocks using options X and written to file
// './target_corrs.dat' (Linux instructions) or './inout/target_corrs.dat'
// (Docker instructions), unless this default name is changed with option C;
// if omitted, target correlations are read from file './target_corrs.dat'
// (Linux instructions) or './inout/target_corrs.dat' (Docker instructions),
// unless this default name is changed with option C
//
// C <string>; name of file where target correlations are read if K is omitted
// and where generated target correlations are written if K is included;
// default: './target_corrs.dat' (Linux instructions) or
// './inout/target_corrs.dat' (Docker instructions)
//
// X<char><char> <double>; identifies target correlations within or between
// subsets of parameters that are randomized (i.e, that vary between
// scenarios); state once for each pair of subsets; target self-correlations
// are equal to 1.0; any use of X without G will be rejected;
// default: zero
// D -> demands
// A -> total capacities of arcs
// C -> commodity-specific capacities of arcs
// F -> fixed costs
// V -> variable costs
// Example: in the instruction line, -XDD 0.9 -XFD 0.1 -XFF 0.01 specifies
// target correlations between distinct demand parameters of 0.9
// (self-correlations being equal to 1.0), between fixed cost parameters and
// demand parameters of 0.1 and between distinct fixed cost parameters of 0.01
// (self-correlations being equal to 1.0); Remark: -XFD and -XDF are equivalent
//
// G <> (no value supplied); if included, target moments are generated based on
// option T and written to file './target_moms.dat' (Linux instructions) or
// './inout/target_moms.dat' (Docker instructions), unless this default name is
// changed with option M; if omitted, target moments are read from file
// './target_moms.dat' (Linux instructions) or './inout/target_moms.dat'
// (Docker instructions) unless this default name is changed with option M
//
// M <string>; name of file where target moments are read if option G is
// omitted or where generated target moments are written if option G is
// included; default: './target_moms.dat' (Linux
// instructions) or './inout/target_moms.dat' (Docker instructions)
//
// T <char>; distributional characteristics of generated target moments
// (U : UNIFORM , T : TRIANGULAR); any use of option T without option G
// will be rejected; default: U
//
// A <double>; parameter alpha used in generating target moments (see
// discussion immediately below); without option G, any use of option A will be
// ignored; default: 0.25
//
// B <double>; parameter beta used in generating target moments (see
// discussion immediately below); without option G, any use of option B will be
// ignored; default: 0.25
//
// When option T is set to U, then the stochastic elements of the network that
// are changed according to scenarios are each assumed to follow a
// uniform(a, b) distribution, where
// a = D - (alpha * D),
// b = D + (beta * D),
// D is the base value of the stochastic element taken from the base
// deterministic network,
// alpha must be in [0, 0.99],
// beta >= 0.
// The distribution is symmetric around D when alpha = beta.
//
// When option T is set to T, then the stochastic elements of the network that
// are changed according to scenarios are each assumed to follow a
// triangular(a, b, c) distribution, where
// a = D - (alpha * D),
// b = D + (beta * D),
// c = D,
// D is the base value of the stochastic element taken from the base
// deterministic network,
// alpha must be in [0, 0.99],
// beta >= 0.
// The distribution is symmetric around D when alpha = beta.
//
// 7.2.2 Options supplied to HKW algorithm
// ***************************************
//
// Scenario generation code is based on paper 'A Heuristic for Moment-
// matching Scenario Generation' by K. Høyland, M. Kaut & S.W. Wallace;
// Computational Optimization and Applications, 24 (2-3), pp. 169–185, 2003;
// doi:10.1023/A:1021853807313. Code by Michal Kaut (michal.kaut@iot.ntnu.no)
// & Diego Mathieu. Additional details about HKW options available in
// ../HOYLAND_KAUT_WALLACE/sg_HKW.c.
//
// n <int>; number of scenarios to generate; default: 1000
//
// t <int>; maximum number of trials (attempts to generate scenarios using
// alternative random starting values); default: 20
//
// i <int>; maximum number of iterations in a trial; default: 50
//
// m <double>; maximum error in matching moments (scaled to var=1);
// default: 0.001
//
// c <double>; maximum error in matching correlations (scaled to var=1);
// default: 0.001
//
// l <int>; level of on-screen reporting by HKW algorithm (between 0 and 11);
// default: 3
//
// o <string>; name of file where resulting matrix of raw scenarios in HKW
// format is written; default: './raw_scens.dat' (Linux instructions) or
// './inout/raw_scens.dat' (Docker instructions)
//
// s <string>; if included, name of file where a matrix of starting scenarios
// saved in HKW format in a previous run (possibly using option o) is read and
// used for starting values; no default
//
// p <string>; if included, read probabilities from this file; otherwise,
// scenarios are assumed equiprobable; no default
//
// r <int>; random seed; default: 7654321
//
// R <int>; random stream; default: 1000
//
// 8 Examples
// **********
//
// The program/examples subfolder holds examples of configuration files and 
// instruction lines occurring when the executable of the program is run
// directly from the Linux instruction line. It also contains the output files
// resulting from each pair of configuration file and instruction line.
//
// Similarly, the docker/examples subfolder holds examples of configuration
// files and instruction lines occurring when the executable of the program is
// run through Docker instructions. These examples are identified by the
// addition of the prefix 'docker_'. The subfolder also contains the output
// files resulting from each pair of configuration file and instruction line.
// The former are identical to the files generated when the executable is run
// from the Linux instruction line.
//
// 8.1 One source and one sink for each commodity
// **********************************************
//
// 8.1.1 Folder ./examples/oneSrcOneSnk/DOW/
// *****************************************
//
// In the base network, there is one source and one sink for each commodity and
// no differentiation of variable costs and capacities according to
// commodities. The base network is received in the DOW format and the
// generated networks are saved in the DOW format.
//
// 8.1.1.1 Operation of the stochastic generator
// *********************************************
//
// Linux command line instruction:
//
// ./stoch_gen +v +F stoch_config_oneSrcOneSnk_10oct2025-18h10_dow.par
//
// Configuration file (input) used with Linux command line:
// stoch_config_oneSrcOneSnk_10oct2025-18h10_dow.par
//
// Docker instruction:
//
// docker run --rm --mount type=bind,src=absolutePathToDockerFolder/inout,dst=/inout stoch_gen_image \
// /stoch_gen +v +F docker_stoch_config_oneSrcOneSnk_10oct2025-18h10_dow.par
//
// Configuration file (input) used with Docker instruction:
// docker_stoch_config_oneSrcOneSnk_10oct2025-18h10_dow.par
//
// Base network file (input):
// determ_oneSrcOneSnk_10oct2025-18h10.dow
//
// Generated moments (output):
// stoch_moms_oneSrcOneSnk_10oct2025-18h10_dow.dat
//
// Generated correlations (output):
// stoch_corrs_oneSrcOneSnk_10oct2025-18h10_dow.dat
//
// Raw HKW scenarios (output):
// stoch_rawscen_oneSrcOneSnk_10oct2025-18h10_dow.dat
//
// Stacked scenarios of the MCFND instance in DOW fomat (output):
// stoch_oneSrcOneSnk_10oct2025-18h10_dow.txt
//
// 8.1.1.2 Origin of the base deterministic network
// ************************************************
//
// The base deterministic network was originally produced by the determistic
// generator as follows.
//
// Linux command line instruction:
//
// ./determ_gen +v +F determ_config_oneSrcOneSnk_10oct2025-18h10.par
//
// Configuration file:
// determ_config_oneSrcOneSnk_10oct2025-18h10.par
//
// Outputs in DOW, STD, LP, MPS, GRA (basic graph), TEX (code of graphical
// representation), pdf (compiled image) formats:
// determ_oneSrcOneSnk_10oct2025-18h10.dow
// determ_oneSrcOneSnk_10oct2025-18h10.std
// determ_oneSrcOneSnk_10oct2025-18h10.mps
// determ_oneSrcOneSnk_10oct2025-18h10.lp
// determ_oneSrcOneSnk_10oct2025-18h10.gra
// determ_oneSrcOneSnk_10oct2025-18h10.tex
// determ_oneSrcOneSnk_10oct2025-18h10.pdf
//
// 8.1.2 Folder ./examples/oneSrcOneSnk/RSTD/
// ******************************************
//
// Similar to 8.1.1, except that the base network is received in the STD format
// and the generated networks are saved in the RSTD format.
//
// 8.1.2.1 Operation of the stochastic generator
// *********************************************
//
// Linux command line instruction:
//
// ./stoch_gen +v +F stoch_config_oneSrcOneSnk_10oct2025-18h10_rstd.par
//
// Configuration file (input) used with Linux command line:
// stoch_config_oneSrcOneSnk_10oct2025-18h10_rstd.par
//
// Docker instruction:
//
// docker run --rm --mount type=bind,src=absolutePathToDockerFolder/inout,dst=/inout stoch_gen_image \
// /stoch_gen +v +F docker_stoch_config_oneSrcOneSnk_10oct2025-18h10_rstd.par
//
// Configuration file (input) used with Docker instruction:
// docker_stoch_config_oneSrcOneSnk_10oct2025-18h10_rstd.par
//
// Base network file (input):
// determ_oneSrcOneSnk_10oct2025-18h10.std
//
// Generated moments (output):
// stoch_moms_oneSrcOneSnk_10oct2025-18h10_rstd.dat
//
// Generated correlations (output):
// stoch_corrs_oneSrcOneSnk_10oct2025-18h10_rstd.dat
//
// Raw HKW scenarios (output):
// stoch_rawscen_oneSrcOneSnk_10oct2025-18h10_rstd.dat
//
// Stacked scenarios of the MCFND instance in RSTD fomat (output):
// stoch_oneSrcOneSnk_10oct2025-18h10_rstd.txt
//
// 8.1.2.2 Origin the base deterministic network
// *********************************************
//
// Similar to 8.1.1.2.
//
// 8.1.3 Folder ./examples/oneSrcOneSnk/GSTD/
// ******************************************
//
// Similar to 8.1.1, except that the base network is received in the STD format
// and the generated networks are saved in the GSTD format (which is identical
// to RSTD when there are one source and one sink for each commodity).
//
// 8.1.3.1 Operation of the stochastic generator
// *********************************************
//
// Linux command line instruction:
//
// ./stoch_gen +v +F stoch_config_oneSrcOneSnk_10oct2025-18h10_gstd.par
//
// Configuration file (input):
// stoch_config_oneSrcOneSnk_10oct2025-18h10_gstd.par
//
// Docker instruction:
//
// docker run --rm --mount type=bind,src=absolutePathToDockerFolder/inout,dst=/inout stoch_gen_image \
// /stoch_gen +v +F docker_stoch_config_oneSrcOneSnk_10oct2025-18h10_gstd.par
//
// Configuration file (input) used with Docker instruction:
// docker_stoch_config_oneSrcOneSnk_10oct2025-18h10_gstd.par
//
// Configuration file (input) used with Linux command line:
// determ_oneSrcOneSnk_10oct2025-18h10.std
//
// Generated moments (output):
// stoch_moms_oneSrcOneSnk_10oct2025-18h10_gstd.dat
//
// Generated correlations (output):
// stoch_corrs_oneSrcOneSnk_10oct2025-18h10_gstd.dat
//
// Raw HKW scenarios (output):
// stoch_rawscen_oneSrcOneSnk_10oct2025-18h10_gstd.dat
//
// Stacked scenarios of the MCFND instance in GSTD fomat (output):
// stoch_oneSrcOneSnk_10oct2025-18h10_gstd.txt
//
// 8.1.3.2 Origin the base deterministic network
// *********************************************
//
// Similar to 8.1.1.2.
//
// 8.2 Multiple but equal sources and multiple but equal sinks
// ***********************************************************
//
// 8.2.1 Folder ./examples/eqSrcsEqSnks/GSTD/
// ******************************************
//
// In the base network, there are multiple sources and multiple sinks for each
// commodity. However, sources and sinks are identical across commodities.
// There is differentiation of variable costs and capacities according to
// commodities. The base network must be received in the STD format (as the DOW
// format is unsuitable in this situation) and the generated networks must be
// saved in the GSTD format.
//
// 8.2.1.1 Operation of the stochastic generator
// *********************************************
//
// Linux command line instruction:
//
// ./stoch_gen +v +F stoch_config_eqSrcsEqSnks_10oct2025-18h10_gstd.par
//
// Configuration file (input) used with Linux command line:
// stoch_config_eqSrcsEqSnks_10oct2025-18h10_gstd.par
//
// Docker instruction:
//
// docker run --rm --mount type=bind,src=absolutePathToDockerFolder/inout,dst=/inout stoch_gen_image \
// /stoch_gen +v +F docker_stoch_config_eqSrcsEqSnks_10oct2025-18h10_gstd.par
//
// Configuration file (input) used with Docker instruction:
// docker_stoch_config_eqSrcsEqSnks_10oct2025-18h10_gstd.par
//
// Base network file (input):
// determ_eqSrcsEqSnks_10oct2025-18h10.std
//
// Generated moments (output):
// stoch_moms_eqSrcsEqSnks_10oct2025-18h10_gstd.dat
//
// Generated correlations (output):
// stoch_corrs_eqSrcsEqSnks_10oct2025-18h10_gstd.dat
//
// Raw HKW scenarios (output):
// stoch_rawscen_eqSrcsEqSnks_10oct2025-18h10_gstd.dat
//
// Stacked scenarios of the MCFND instance in DOW fomat (output):
// stoch_eqSrcsEqSnks_10oct2025-18h10_gstd.txt
//
// 8.2.1.2 Origin of the base deterministic network
// ************************************************
//
// The base deterministic network was originally produced by the determistic
// generator as follows.
//
// Linux command line instruction:
//
// ./determ_gen +v +F determ_config_eqSrcsEqSnks_10oct2025-18h10.par
//
// Configuration file:
// determ_config_eqSrcsEqSnks_10oct2025-18h10.par
//
// Outputs in DOW, STD, LP, MPS, GRA (basic graph), TEX (code of graphical
// representation), pdf (compiled image) formats:
// determ_eqSrcsEqSnks_10oct2025-18h10.dow
// determ_eqSrcsEqSnks_10oct2025-18h10.std
// determ_eqSrcsEqSnks_10oct2025-18h10.mps
// determ_eqSrcsEqSnks_10oct2025-18h10.lp
// determ_eqSrcsEqSnks_10oct2025-18h10.gra
// determ_eqSrcsEqSnks_10oct2025-18h10.tex
// determ_eqSrcsEqSnks_10oct2025-18h10.pdf
//
//
#include <HOYLAND_KAUT_WALLACE/matrix.h>
#include <configuration.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

#include "./instance.h"
#include "./model.h"

using namespace std;

typedef vector<int> VectInt;

const int numbMoments = 4;

int CommodOnArc::count = 0;
double *CommodOnArc::ptrCap = NULL;
double *CommodOnArc::ptrCost = NULL;

int Commodity::count = 0;
double *Commodity::ptr = NULL;

int Arc::count = 0;
double *Arc::ptrCap = NULL;
double *Arc::ptrCost = NULL;

int main(int argc, char **argv)
{

    // process configuration options
    Configuration params(argc, argv, "./stoch_config.par");
    Instance *inst;

    // initialize computational instance
    if (params.instanceFormat == DOW)
        inst = new InstanceDow(&params);
    else if (params.instanceFormat == RSTD)
        inst = new InstanceRStd(&params);
    else if (params.instanceFormat == GSTD)
        inst = new InstanceGStd(&params);

    // create input and output streams
    ifstream is(params.inputInstFileName);

    if (!is.is_open())
    {
        cerr << "\nERROR: " << endl;
        cerr << "The file " << params.inputInstFileName
             << " containing the base deterministic network could not be opened." << endl;
        cerr << "EXECUTION ABORTED\n" << endl;
        exit(1);
    }

    ofstream os(params.outputInstFileName);
	
	if (!os.is_open())
    {
        cerr << "\nERROR: " << endl;
        cerr << "The file " << params.outputInstFileName
             << " where network instances resulting from scenario generation\n"
			 << "should be written could not be opened." << endl;
        cerr << "EXECUTION ABORTED\n" << endl;
        exit(1);
    }

    // read base deterministic network from source;
    // input operator is defined in instance.h;
    // InstanceDow::read() or InstanceStd::read() is called
    is >> *inst;

    // generate all raw scenarios and store in computational instance
    int numbScenarios = inst->getScenarios();

    // create solution model for feasibility testing
    Model *model;

#if CPLEX

    // case GSTD: not one source node and one sink node per commodity
    if (params.instanceFormat == GSTD)
        model = new CPXModel_GStd(inst);
    // cases DOW or RSTD: one source node and one sink node per commodity
    else
        model = new CPXModel_DowRStd(inst);

#elif HIGHS

    // case GSTD: not one source node and one sink node per commodity
    if (params.instanceFormat == GSTD)
        model = new HIGHSModel_GStd(inst);
    // cases DOW or RSTD: one source node and one sink node per commodity
    else
        model = new HIGHSModel_DowRStd(inst);

#endif

    // serves to count feasible scenarios
    int numbFeasible = 0;

    cout << "\nVerifying feasibility of every generated scenario.\n" << endl;

    // iterate over generated raw scenarios
    for (int j = 0; j < numbScenarios; j++)
    {

        if (j % 200 == 0)
            cout << "Currently verifying feasibility of " << j + 1 << "th scenario out of " << numbScenarios << "."
                 << endl;

        // update instance with jth raw scenario
        inst->copyScenario(j);

        // check if instance updated with scenario is feasible when all arcs
        // are opened; write instance to output stream if so
        if (model->testScenario())
        {
            // output operator is defined in instance.h;
            // InstanceDow::write() or InstanceStd::write() is called
            os << *inst << endl;
            numbFeasible++;
        }
    }

    cout << "\nSaved " << numbFeasible << " feasible scenarios out of " << numbScenarios << " generated. " << endl;
    os.close();
    delete model;
    delete inst;
}
