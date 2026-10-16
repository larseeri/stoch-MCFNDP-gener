(October 16th 2026)

## Generator of instances of the linear, two-stage, stochastic, multi-commodity, capacitated, fixed charge network design problem with recourse

## Tutorial example


To become proficient with the use of the GenStoch stochastic generator, we suggest to start with steps 1-2 and to iterate on steps 3-5:

1 Read Sections 1 and 3 of Larsen et al. (2023) attentively.\
2 Skim the user guide USERGUIDE_STOCH.md entirely to build a cognitive structure (know where to seek information, start locating the information related to the task at hand). The user guide details installation, operation, functionalities and available options of the generator.\
3 Study this tutorial example in detail.\
4 Study USERGUIDE_STOCH.md in detail.\
5 Launch generation experiments.


Using the generator involves the following steps A, B and C.

### A- If not already done, install the program or create a Docker container from the supplied Docker image.

See Sections 1 and 2 of the online user guide for complete details.


### B- Build or modify a configuration file.

As an example, we reproduce this configuration file taken from the stochastic library:

library/hubsAndSats/eqSrcsEqSnks/flat/hs_eqSrcsEqSnks_flat_base-det-seed-50163844_base-det-strm-5343_smlInstLoo_F-G_G_T-T_A-0~5_B-0~5_S-3_K_XDD-0~5_XDA-0~2_XAA-0~5_n-3000_t-20_i-50_l-3_m-0~001_c-0~001_seed-50163844_strm-22912_gstd.par

By thoroughly understanding this example and experimenting with alternative and/or additional settings, the user will gain sufficient knowledge to use the full extent of the available functionalities in new experiments.


#### 1 Content of the configuration file

We reordered the content of the original file to make the exposition flow more progressively.

*** Content of the configuration file starts on the next line. ***

// this configuration file is called with this command line:\
// ./stoch_gen +v +F ../library/hubsAndSats/eqSrcsEqSnks/flat/smlInstLoo/hs_eqSrcsEqSnks_flat_base-det-seed-50163844_base-det-strm-5343_smlInstLoo_F-G_G_T-T_A-0~5_B-0~5_S-3_K_XDD-0~5_XDA-0~2_XAA-0~5_n-3000_t-20_i-50_l-3_m-0~001_c-0~001_seed-50163844_strm-22912_gstd.par

// base deterministic instance file has GSTD format\
-F G

// path of base deterministic instance\
-I ../library/hubsAndSats/eqSrcsEqSnks/flat/hs_vols-20-100_csts-20-100-20-100-20-100-20-100_eqSrcsEqSnks_secTop-flat_secTopSpe-2-3-2-3-2-3-2-8-2-8-12-15-45-10-15-2-15-45-10-15-2_noPar-20_adjSpe-0~33-1~00-0~33-1~00_seed-50163844_strm-5343.std

// moments will be generated\
G

// distribution\
-T T

// alpha\
-A 0.5

// beta\
-B 0.5

// randomizing volumes and overall capacities\
-S 3\
// correlations are hereafter specified for each block\
-K\
-XDD 0.5\
-XDA -0.2\
-XAA 0.5

// number of scenarios\
-n 3000

// maximum number of trials in HKW algo.\
-t 20

// maximum number of iterations per trial in HKW algo.\
-i 50

// level of onscreen reporting\
-l 3

// tolerance about errors in moments\
-m 0.001

// tolerance about errors in correlations\
-c 0.001

// seed\
-r 50163844

// stream\
-R 22912

// path of file where correlations will be saved\
-C ../library/hubsAndSats/eqSrcsEqSnks/flat/smlInstLoo/corrs_hs_eqSrcsEqSnks_flat_base-det-seed-50163844_base-det-strm-5343_smlInstLoo_F-G_G_T-T_A-0~5_B-0~5_S-3_K_XDD-0~5_XDA-0~2_XAA-0~5_n-3000_t-20_i-50_l-3_m-0~001_c-0~001_seed-50163844_strm-22912_gstd.dat

// path of file where moments will be saved\
-M ../library/hubsAndSats/eqSrcsEqSnks/flat/smlInstLoo/moms_hs_eqSrcsEqSnks_flat_base-det-seed-50163844_base-det-strm-5343_smlInstLoo_F-G_G_T-T_A-0~5_B-0~5_S-3_K_XDD-0~5_XDA-0~2_XAA-0~5_n-3000_t-20_i-50_l-3_m-0~001_c-0~001_seed-50163844_strm-22912_gstd.dat

// path of file where raw scenarios will be saved\
-o ../library/hubsAndSats/eqSrcsEqSnks/flat/smlInstLoo/raw-scens_hs_eqSrcsEqSnks_flat_base-det-seed-50163844_base-det-strm-5343_smlInstLoo_F-G_G_T-T_A-0~5_B-0~5_S-3_K_XDD-0~5_XDA-0~2_XAA-0~5_n-3000_t-20_i-50_l-3_m-0~001_c-0~001_seed-50163844_strm-22912_gstd.dat

// path of file where generated stochastic instance will be saved\
-O ../library/hubsAndSats/eqSrcsEqSnks/flat/smlInstLoo/stoch-inst_hs_eqSrcsEqSnks_flat_base-det-seed-50163844_base-det-strm-5343_smlInstLoo_F-G_G_T-T_A-0~5_B-0~5_S-3_K_XDD-0~5_XDA-0~2_XAA-0~5_n-3000_t-20_i-50_l-3_m-0~001_c-0~001_seed-50163844_strm-22912_gstd.txt

*** Content of the configuration file ends at the preceding line. ***


#### 2 Detailed explanations

We explain each statement and comment appearing in the example.

Remarks:

a- To supply an option, a line of a configuration file must begin with a dash '-' or a letter. A line beginning with any other character will be considered a comment. It is advisable to start comment lines exclusively with "//" or with "#". Blank lines are ignored. See Sections 1.2, 2 and 7.2 of the online user guide for complete details.

b- The order of the statements and comments in a configuration file is immaterial to the workings of the generator.

c- Any boolean option that does not appear in the configuration file or in the command line is interpreted as having the boolean value false. To set a boolean option to true, it must appear without any value in the configuration file or in the command line. Option values cannot be set explicitly with boolean values true or false. 

d- The selected option values must conform to the restrictions stated in Section 4.3 of the online user guide. The program will notify the user of any error and stop.

e- In this section, the suffixes appearing in the names of the option values indicate the expected type. For instance, option F expects a character. Hence the value is named F_char. 

**2.1** It is advisable to make sure that the name of the configuration file is a unique identifier of its content and to state the command line calling the current configuration file as a comment within the configuration file itself. Here, name_string==./stoch_gen +v +F ../library/hubsAndSats/eqSrcsEqSnks/flat/smlInstLoo/hs_eqSrcsEqSnks_flat_base-det-seed-50163844_base-det-strm-5343_smlInstLoo_F-G_G_T-T_A-0~5_B-0~5_S-3_K_XDD-0~5_XDA-0~2_XAA-0~5_n-3000_t-20_i-50_l-3_m-0~001_c-0~001_seed-50163844_strm-22912_gstd.par.

```
// file name: name_string
```

**2.2** Specifies the format of the base deterministic instance file. Here, F_char=='G' whence it has GSTD format.

```
-F F_char
```

**2.3** Specifies the path of base deterministic instance. See section 5.1 of online user guide. Here,  I_string==../library/hubsAndSats/eqSrcsEqSnks/flat/hs_vols-20-100_csts-20-100-20-100-20-100-20-100_eqSrcsEqSnks_secTop-flat_secTopSpe-2-3-2-3-2-3-2-8-2-8-12-15-45-10-15-2-15-45-10-15-2_noPar-20_adjSpe-0~33-1~00-0~33-1~00_seed-50163844_strm-5343.std. This is the instance described in the commented example available online about the deterministic generator.

```
-I I_string
```

**2.4** Signals that target moments will be generated under a distributional assumption. If absent, target moments will be read from file specified with option M. If present, target moments will be written to the specified file. Here, target moments will be generated. See Section 3.2 of the online user guide.

```
-G
```

**2.5** Specifies which assumed distribution will be used to generate target moments. See Section 3.2 of the online user guide. Here, T_char=='T' whence the assumed distribution is triangular.

```
-T T_char
```

**2.6** Specifies the proportional distance from any element of the base network that will be randomized to the lower bound of its assumed distribution. See Section 3.2 of the online user guide. Here, A_double==0.5, whence the lower bound will be at a value 50% below the base value.

```
-A A_double
```

**2.7** Specifies the proportional distance from any element of the base network that will be randomized to the upper bound of its assumed distribution. See Section 3.2 of the online user guide. Here, B_double==0.5, whence the upper bound will be at a value 50% above the base value.

```
-B B_double
```

**2.8** Specifies the integral values r_int and R_int that are respectively assigned to the seed and stream of the pseudo-random generator. To ensure reproducibility, they should always be set. When several configurations are included in a set, the pairs (seed_number, stream_number) should be selected for each configuration so as to constitute a unique identifier. A practical choice is to make seed_number invariant within the set of configurations whereas stream_number varies among the configuration in an interpretable way. Here, r_int==50163844 and R_int==2291**2.

```
-r r_int
```

```
-R r_int
```

**2.9** Specifies which elements of the of the base deterministic network must be randomized. See Section 3.1 of the online user guide. Here, S_int==3 and since 3=1+2, demands and overall capacities are to be randomized.

```
-S s_int
```

**2.10** Signals that blocks of the target correlation matrix of randomized elements will be assigned values. If absent, the target correlation matrix will be read from the file specified with option C. If present, target correlations will be written to specified file. See Section 3.3 of the online user guide. Here, blocks of the target correlation matrix will be assigned values.

```
-K
```

**2.11** Specify the values assigned to blocks of the target correlation matrix. These blocks are identified from option S. See Section 3.3 of the online user guide. XDD XDD_double indicates that target cross-correlations among demands will be equal to XDD_double; XDA XDA_double indicates that target cross-correlations between demands and overall capacities will be equal to XDA_double; XAA XAA_double indicates that target cross-correlations among overall capacities will be equal to XAA_double. Here, XDD_double==0.5, XDA_double==-0.2, XAA_double==0.5.

```
-XDD XDD_double
-XDA XDA_double
-XAA XAA_double
```

**2.12** Specifies the number of scenarios. Here, n_int==3000.

```
-n n_int
```

**2.13** (Not used in the current example) Specifies the path where the probabilities are to be read. If absent, equiprobability is assumed.

```
-p p_string
```

**2.14** Specifies the maximum number of trials in HKW algo. Here, t_int==20.

```
-t t_int
```

**2.15** Specifies the maximum number of iterations per trial in HKW algo. Here, i_int==50.

```
-i i_int
```

**2.16** Specifies the level of onscreen reporting. Here, l_int==3.

```
-l l_int
```

**2.17** Specifies the tolerance about RMS error between target and realized standardized moments. Here, m_double==0.001.

```
-m m_double
```

**2.18** Specifies the tolerance about RMS error between target and realized correlations. Here, c_double==0.001.

```
-c c_double
```

**2.19** Specifies path of file where target correlations will be saved or read. See Section 5.2.2 of online user guide.v Here, target correlations are written to C_string==../library/hubsAndSats/eqSrcsEqSnks/flat/smlInstLoo/corrs_hs_eqSrcsEqSnks_flat_base-det-seed-50163844_base-det-strm-5343_smlInstLoo_F-G_G_T-T_A-0~5_B-0~5_S-3_K_XDD-0~5_XDA-0~2_XAA-0~5_n-3000_t-20_i-50_l-3_m-0~001_c-0~001_seed-50163844_strm-22912_gstd.dat.

```
-C C_string 
```

**2.20** Specifies path of file where target moments will be saved or read. See Section 5.2.1 of online user guide. Here, target moments are written to M_string==../library/hubsAndSats/eqSrcsEqSnks/flat/smlInstLoo/moms_hs_eqSrcsEqSnks_flat_base-det-seed-50163844_base-det-strm-5343_smlInstLoo_F-G_G_T-T_A-0~5_B-0~5_S-3_K_XDD-0~5_XDA-0~2_XAA-0~5_n-3000_t-20_i-50_l-3_m-0~001_c-0~001_seed-50163844_strm-22912_gstd.dat.

```
-M M_string
```

**2.21** Specifies path of file where raw scenarios will be saved. See Section 5.2.4 of online user guide. Here, o_string==../library/hubsAndSats/eqSrcsEqSnks/flat/smlInstLoo/raw-scens_hs_eqSrcsEqSnks_flat_base-det-seed-50163844_base-det-strm-5343_smlInstLoo_F-G_G_T-T_A-0~5_B-0~5_S-3_K_XDD-0~5_XDA-0~2_XAA-0~5_n-3000_t-20_i-50_l-3_m-0~001_c-0~001_seed-50163844_strm-22912_gstd.dat.

```
-o o_string 
```

**2.22** Specifies path of file where the generated stochastic instance will be saved. See Section 4 of online user guide. Here, O_string==../library/hubsAndSats/eqSrcsEqSnks/flat/smlInstLoo/stoch-inst_hs_eqSrcsEqSnks_flat_base-det-seed-50163844_base-det-strm-5343_smlInstLoo_F-G_G_T-T_A-0~5_B-0~5_S-3_K_XDD-0~5_XDA-0~2_XAA-0~5_n-3000_t-20_i-50_l-3_m-0~001_c-0~001_seed-50163844_strm-22912_gstd.txt.

```
-O O_string
```

**2.23** (Unused in the current example) Signals that all scenarios of the resulting stochastic MCFNDP instance are to be written to file without first checking for feasibility of each scenario and reporting on results. (Useful when regenerating an instance of a library whose scenarios' feasibility has already been demonstrated.)

```
-J
```


### C- Starting the generation process

See Sections 1-2 of the online user guide for a detailed discussion.

For example, to launch the generation process from the Linux command line based on the configuration file described in Section B above, the following instruction would be issued from the folder /program.

```
./stoch_gen +v +F ../library/hubsAndSats/eqSrcsEqSnks/flat/smlInstLoo/hs_eqSrcsEqSnks_flat_base-det-seed-50163844_base-det-strm-5343_smlInstLoo_F-G_G_T-T_A-0~5_B-0~5_S-3_K_XDD-0~5_XDA-0~2_XAA-0~5_n-3000_t-20_i-50_l-3_m-0~001_c-0~001_seed-50163844_strm-22912_gstd.par
```

Notice that options can also be set on the command line. Furthermore, the value of an option specified on the command line supersedes that specified in the configuration file.
