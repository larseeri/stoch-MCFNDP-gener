EXT = cpp
CXX = g++
EXEC = stoch_gen

CPLEX_SOLVER = 0
HIGHS_SOLVER = 0

# CPLEX locations
CPLEX_HOME   = /home/ibm/cplex-studio/22.1.1/
CPLEX_INCDIR = $(CPLEX_HOME)/cplex/include
CONCERT_INCDIR = $(CPLEX_HOME)/concert/include
CPLEX_LIBDIR   = $(CPLEX_HOME)/cplex/lib/x86-64_linux/static_pic/
CONCERT_LIBDIR   = $(CPLEX_HOME)/concert/lib/x86-64_linux/static_pic/

# HIGHS locations
HIGHS_HOME   = ./HIGHS/HiGHS-1.11.0/install
HIGHS_INCDIR = $(HIGHS_HOME)/include/highs/
HIGHS_LIBDIR   = $(HIGHS_HOME)/lib/

# internal locations
HOYLAND_DIR  = ./HOYLAND_KAUT_WALLACE
RANDOM_DIR   = ./PCG_RANDOM/include
GENSTOCH_DIR = ./GENSTOCH

ADDFLAGS = -g -Wall #-Wextra #-O3

ifeq ($(CPLEX_SOLVER), 1)
    CXXFLAGS = -DIL_STD -DNDEBUG -DCPLEX
    INCFLAGS = -I. -I${HOYLAND_DIR} -I${RANDOM_DIR} -I${GENSTOCH_DIR} -I$(CPLEX_INCDIR) -I$(CONCERT_INCDIR)
    LDFLAGS  = -L. -L$(CPLEX_LIBDIR) -L $(CONCERT_LIBDIR) -lilocplex -lcplex -lconcert -lm -lpthread -ldl
else 
    ifeq ($(HIGHS_SOLVER), 1)
        CXXFLAGS = -DHIGHS
        INCFLAGS = -I. -I${HOYLAND_DIR} -I${RANDOM_DIR} -I${GENSTOCH_DIR} -I$(HIGHS_INCDIR)
        LDFLAGS  = -L. -L$(HIGHS_LIBDIR) -lhighs -Wl,-rpath,$(HIGHS_LIBDIR)
    else
        $(error ERROR: Specify either CPLEX_SOLVER=1 or HIGHS_SOLVER=1)
    endif
endif

SRC = $(wildcard ${GENSTOCH_DIR}/*.cpp)
SRC2 = $(wildcard ${HOYLAND_DIR}/*.c)
OBJ = $(SRC:.cpp=.o)
OBJ2 = $(SRC2:.c=.o)
OBJDIR = ${GENSTOCH_DIR}
OBJDIR2 = ${HOYLAND_DIR}

all: $(EXEC)

$(EXEC): $(OBJ) $(OBJ2)
	g++ -o $@ $^  $(OPTFLAGS) $(LDFLAGS) $(ADDFLAGS)

$(OBJDIR)/%.o: ${GENSTOCH_DIR}/%.cpp ${GENSTOCH_DIR}/*.h
	g++ -o $@ -c $<  $(ADDFLAGS) $(CXXFLAGS) $(INCFLAGS) $(LDFLAGS)

$(OBJDIR2)/%.o: ${HOYLAND_DIR}/%.c ${HOYLAND_DIR}/*.h
	g++ -o $@ -c $<  $(ADDFLAGS) $(CXXFLAGS) $(INCFLAGS) $(LDFLAGS)

clean:
	@rm -rf *.o ${HOYLAND_DIR}/*.o ${GENSTOCH_DIR}/*.o
	@rm -f $(EXEC)
