.SUFFIXES: .F90 .o
TARGET = mt_cafe


# fortran90/95 compiler
FC  := ifort
#FC := gfortran

# fortran compiler options
FLAGS := -O3

MTMOD = mt_kind_defs.o mt_stream.o gf2xe.o
MTSUB = f_get_coeff.o


all : $(TARGET)

$(TARGET) : $(MTMOD) $(MTSUB) mt_cafe.o
	$(FC) $(FLAGS) $^ -o $@

clean :
	rm -f *.o *.mod $(TARGET)

# suffix rule
.F90.o:
	$(FC) $(FLAG) -c $< 
