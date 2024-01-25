Copyright (c) 2023 Eric Larsen

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the “Software”), to
deal in the Software without restriction, including without limitation the
rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
sell copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
IN THE SOFTWARE.
****************************************************************************
****************************************************************************

### Scenarios generator for the linear, two-stage, multi-<br/>commodity, fixed charge, capacitated network design<br/>problem with stochastic recourse: User's guide (19Nov2023)

This application randomly generates scenarios by matching target moments and
correlations for specified subsets of parameters of a given instance of the linear
deterministic multi-commodity, fixed charge, capacitated network design (MCFND) problem. Details of this process are governed by user-specified configuration settings. The program published by Kaut that implements the algorithm of Hoyland, Kaut and Wallace (HKW) is used to compute the scenarios.

For detailed background, see
Larsen, E., Bisaillon, S., Cordeau, J.F. and Frejinger, E. 
Pseudo-random Instance Generators in C++ for Deterministic and Stochastic 
Multi-commodity Network Design Problems.  

These two paragraphs are adapted from ../HOYLAND_KAUT_WALLACE/sg_HKW.c:

Scenario generation code is based on paper 'A Heuristic for Moment-matching
Scenario Generation' by K. Høyland, M. Kaut & S.W. Wallace; Computational
Optimization and Applications, 24 (2-3), pp. 169–185, 2003;
doi:10.1023/A:1021853807313. Code by Michal Kaut (michal.kaut@iot.ntnu.no)
& Diego Mathieu.

The HKW code generates scenarios for multivariate random variables by
matching specified target first four moments and correlations. The algorithm
is iterative and reports distance from target values. Convergence is not
guaranteed and the algorithm can run in several trials if needed, starting
each one from new random initializations.

#### 1 Building the program

Under Linux, follow these steps to build the program:

a- Install CPLEX (if not already done) and note the installation directory.<br/>
b- Uncompress the archive file containing the program files.<br/>
c- Go to the directory just created.<br/>
d- Edit file 'makefile' to set the variable CPLEX_HOME to the one noted in
	step a- and save modified makefile to disk.<br/>
e- In the terminal, run this command:  make
	to compile the program.<br/>
f- An executable file named 'exe' will be created.

#### 2 Running the program

The general command for running the program on Linux is as follows:

./exe -S <stochasticElements> [main options] [HKW options]

If running ./exe without any options, the list of all options for the
managing code (identified with capital letters) and for the HKW program
(identified with lower case letters) will be displayed.

Option '-S <number>' must be invoked to determine which subsets of parameters
are meant to vary between scenarios (we call them the stochastic
elements). This option (-S) expects a number (from 1 to 31) that is a
sum of the following bits:<br/>
1 -> demands<br/>
2 -> total capacities of arcs<br/>
4 -> commodity-specific capacities of arcs<br/>
8 -> fixed costs<br/>
16 -> variable costs

Other options have default values that are hard-coded in the file
configuration.cpp.

Hence, running this command line:

./exe -S 11

indicates to the program that fixed costs (8), total capacities of arcs (2) and
demands (1) must be randomized (8 + 2 + 1 = 11).

The value of any option specified on the command line will override that
specified in file configuration.cpp. If the value of an option is specified
more than once on the command line, the last value is retained.

DETAILED USAGE OF ALL OPTIONS IS SUPPLIED FURTHER BELOW.

Here are some pointers and two examples:

The default input instance is set to test.dow and the default input format is
set to DOW. If the input instance is changed with option
'-I <input instance>', it is critical that the input format is set
accordingly with option '-F <format>'. For example:

./exe -S 3 -I instB.std -F S

indicates that demands and total capacities of arcs will be randomized based
on the instance instB.std which has the format 'S' (thus STD).

Target moments and correlations can be EITHER specified in files OR
generated by the program based on default settings and indications supplied
on the command line. If they are to be generated, use option '-G'.
A single target correlation coefficient within subset of parameters U or
between subsets of parameters U and V can be specified with '-X<UU> <value>'
and '-X<UV> <corr>' respectively, where corr is the target correlation
coefficient and U and V can take these values:<br/>
D -> demands<br/>
A -> total capacities of arcs<br/>
C -> commodity-specific capacities of arcs<br/>
F -> fixed costs<br/>
V -> variable costs<br/>
Target correlations that are unspecified are assumed by default to be zero.
Target moments are generated based (i) on a distributional assumption that
can be specified using option '-T <dist>' where dist stands for either U
(uniform) or D (triangular) and (ii) on two parameters, alpha and beta, that
can be specified using options '-A <alpha>' and '-B <beta>'. Any use of
options -X, -T, -A or -B without option -G will be rejected. Hence, the
following command line specifies target correlations among demands of 0.5,
between demands and total capacities of arcs of -0.3, and among total
capacities of arcs of 0.7. It specifies that the target first four moments
of all randomized parameters will be calculated based on uniform
distributions using alpha = beta = 0.25.

ex.:  ./exe -S 3 -I instB.std -F S -G -XDD 0.5 -XDA -0.3 -XAA 0.7 -T U -A 0.25 -B 0.25

#### 3 Output file

The output file yielded by the program superposes complete deterministic
problem instances, one for each scenario realization. Each instance is
preceded by a separator as follows:

SCENARIO   <scenario number>

Immediately after the separator, the instance corresponding to the
identified scenario is fully described with the same format, DOW or STD, as
that used to supply the base deterministic problem instance. Stacking and
separating full descriptions of scenario realizations in this fashion makes
it easy to pinpoint and read the first and second stages information
relevant to every individual scenario.

#### 4 Format of input files supplied to HKW algorithm

(adapted from ../HOYLAND_KAUT_WALLACE/sg_HKW.c)

##### 4.1 Target moments and correlations

When option -G is absent from command line and unless overridden by parameter
-M, target moments are read from file named 'tg_moms.txt'.

When option -G is absent from command line and unless overridden by parameter
-C, target correlations are read from file named 'tg_corrs.txt'.

These text files include a matrix of numbers in the following format
(hence the target moment file has 6 rows):<br/>
number of rows<br/>
number of columns<br/>
data (by rows)<br/>
(rest of the file is ignored)

IMPORTANT: Expected definitions of the 4 target moments are as follows:<br/>
(i) population estimators are used (as in spreadsheets)<br/>
(ii) 2nd moment is variance and NOT standard deviation<br/>
(iii) 4th moment is kurtosis - 3<br/>
(iv) 3rd and 4th moments are normalized with standard deviation

##### 4.2 Probabilities of scenarios

If the probabilities ARE NOT supplied in a file, then scenarios are assumed
to be equiprobable. If the probabilities ARE supplied in a file, the file has
to be in a vector format:<br/>
number of elements<br/>
data<br/>
(rest of the file is ignored)

##### 4.3 Scenarios in HKW format

Scenarios saved in HKW format in a previous run (using option -o) can be
used as starting values in a next run (using option -s). These will only
be used in the first trial. If this does not converge and the maximum
number of trials specified by option -t is larger than 1, the next trials
start from a random starting point.<br/>
IMPORTANT: Scenarios in HKW format written in file identified with option -o
are distinct from desired output of the program written in file identified
with option -O and described in Section 3 above.

#### 5 Usage of options

(some parts adapted from ../HOYLAND_KAUT_WALLACE/sg_HKW.c)

##### 5.1 Main generator options

-h <> (no value supplied); display help message

-F <char>; format of base deterministic input instance instance
(D: DOW, S: STD); default: D

DOW and STD formats are specific to MCFND problems and described in files
DOW-format-desc.txt and STD-format-desc.txt; DOW presupposes that for each
arc, capacities and variable costs are identical for all commodities;
LP supplies a generic human readable statement

-I <string>; name of file where deterministic MCFND base instance is to be
read; default: test.dow

-O <string>; name of file where instances resulting from scenario
generations are to be written; default: output.txt

-S <int>; identifies which subsets of parameters are meant to vary between
scenarios; expects a number (from 1 to 31) that is a sum of the following
bits:<br/>
1 -> demands<br/>
2 -> total capacities of arcs<br/>
4 -> commodity-specific capacities of arcs<br/>
8 -> fixed costs<br/>
16 -> variable costs

-G <> (no value supplied); if present, target moments and target
correlations are generated based on options -T and -X and WRITTEN according
to options -M and -C; if absent, they are READ according to options -M and
-C.

-C <string>; name of file where target correlations are read if -G is absent
and where generated target correlations are written if -G is present;
default: tg_corrs.txt

-X<char><char> <double>; identifies target correlations within sets of
parameters or between sets of parameters (state once for each pair of sets):<br/>
D -> demands<br/>
A -> total capacities of arcs<br/>
C -> commodity-specific capacities of arcs<br/>
F -> fixed costs<br/>
V -> variable costs<br/>
example: -XDD 0.9 -XFD 0.1 -XFF 0.01 specifies target correlations within demand parameters of 0.9, between fixed cost parameters and demand parameters
of 0.1 and within fixed cost parameters of 0.01; any use of -X without -G will be rejected; N.B.: -XFD and -XDF are equivalent; default: zero correlation

-M <string>; name of file where target moments are read if -G is absent or
where generated target moments are written if -G is present;
default: tg_moms.txt

-T <char>; distributional characteristics of generated target moments
( U : UNIFORM , T : TRIANGULAR ); any use of -T without -G will be rejected; 
default: U

-A <double>; parameter alpha used in generating target moments (see below); 
any use of -A without -G will be rejected; 
default: 0.25

-B <double>; parameter beta used in generating target moments (see below);<br/>
any use of -B without -G will be rejected; 
default: 0.25

When MOMENTS_DIST_CHARACT == UNIFORM, then the stochastic elements
of the network that are changed according to scenarios are each assumed
to follow a uniform(a, b) distribution where,<br/>
a = D - (alpha * D),<br/>
b = D + (beta * D),<br/>
D is the base value of the stochastic element taken from the base
deterministic network,<br/>
alpha in [0, 1),<br/>
beta >= 0.<br/>
The distribution is symmetric around D when alpha = beta.

When MOMENTS_DIST_CHARACT == TRIANGULAR, then the stochastic elements
of the network that are changed according to scenarios are each assumed
to follow a triangular(a, b, c) distribution where,<br/>
a = D - (alpha * D),<br/>
b = D + (beta * D),<br/>
c = D,<br/>
D is the base value of the stochastic element taken from the base
deterministic network,<br/>
alpha in [0, 1),<br/>
beta >= 0.<br/>
The distribution is symmetric around D when alpha = beta.

##### 5.2 Options supplied to HKW algorithm

Scenario generation code is based on paper 'A Heuristic for Moment-matching
Scenario Generation' by K. Høyland, M. Kaut & S.W. Wallace; Computational
Optimization and Applications, 24 (2-3), pp. 169–185, 2003;
doi:10.1023/A:1021853807313. Code by Michal Kaut (michal.kaut@iot.ntnu.no)
& Diego Mathieu. Additional details about HKW options in
../HOYLAND_KAUT_WALLACE/sg_HKW.c.

-n <int>; number of scenarios to generate; default: 100

-t <int>; maximum number of trials (attempts to generate) using alternative
random starting values; default: 10

-i <int>; maximum number of iterations in a trial; default: 20

-m <double>; maximal error in matching moments (scaled to var=1);
default: 0.001

-c <double>; maximal error in matching correlations (scaled to var=1);
default: 0.001

-l <int>; level of on-screen reporting by HKW algorithm (between 0 and 11);
default: 2

-o <string>; name of file where resulting matrix of scenarios in HKW format
is written; default: out_scen.txt

-s <string>; if present, name of file where a matrix of starting scenarios
saved in HKW format in a previous run (using option -o) is read and used for
starting values

-p <string>; if present, read probabilities from this file; otherwise,
scenarios are assumed equiprobable

-r <int>; random seed; default: 0

-R <int>; random stream; default: 1

Example of command line:
./exe -I test.dow -F D -G -i 500 -S 3 -XDD 0.5 -XAA 0.5 -XAD -0.2 -n 10000 -R 1212121 -r 1212 -A 0.5 -B 0.5
