# programming environment
COMPILER     := mpic++
INCLUDE      := 
LIB          :=  -lfftw3 -lm -lhdf5 -lgsl -lgslcblas -lclass  
# target and source
EXEC         := kevolution
SOURCE       := main.cpp
HEADERS      := $(wildcard *.hpp)

# mandatory compiler settings (LATfield2)
DLATFIELD2   := -DFFT3D -DHDF5

# optional compiler settings (LATfield2)
#DLATFIELD2   += -DH5_HAVE_PARALLEL
#DLATFIELD2   += -DEXTERNAL_IO # enables I/O server (use with care)
#DLATFIELD2   += -DSINGLE      # switches to single precision, use LIB -lfftw3f

# optional compiler settings (gevolution)
DGEVOLUTION  := -DPHINONLINEAR
DGEVOLUTION  += -DBENCHMARK
#DGEVOLUTION  += -NONLINEAR_TEST # for the non-linear instability tests
DGEVOLUTION  += -DEXACT_OUTPUT_REDSHIFTS
#DGEVOLUTION  += -DVELOCITY      # enables velocity field utilities
DGEVOLUTION  += -DCOLORTERMINAL
#DGEVOLUTION  += -DCHECK_B
DGEVOLUTION  += -DHAVE_HICLASS    # -DHAVE_HICLASS  or -DHAVE_CLASS requires LIB -lclass. The initial conditions are provided by hiclass!
DGEVOLUTION  += -DHAVE_HICLASS_BG    # -DHAVE_HICLASS requires LIB -lclass. The BG quantities are provided by hiclass and also parameters like c_s^2,w ...
#DGEVOLUTION  += -DHAVE_HEALPIX  # requires LIB -lchealpix

CDBG +=
CFLAGS += $(CDBG)

# further compiler options
OPT          := -O3 -std=c++11

$(EXEC): $(SOURCE) $(HEADERS) makefile
	$(COMPILER) $< -o $@ $(OPT) $(DLATFIELD2) $(DGEVOLUTION) $(INCLUDE) $(LIB)
	
lccat: lccat.cpp
	$(COMPILER) $< -o $@ $(OPT) $(DGEVOLUTION) $(INCLUDE)
	
lcmap: lcmap.cpp
	$(COMPILER) $< -o $@ $(OPT) -fopenmp $(DGEVOLUTION) $(INCLUDE) $(LIB) $(HPXCXXLIB)

clean:
	-rm -f $(EXEC) lccat lcmap

