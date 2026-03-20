FC = gfortran
FFLAGS = -fopenmp -O3
LDFLAGS = -fopenmp

TARGET = solver

SRC = precision_mod.f90 linear_solver.f90 main.f90
OBJ = $(SRC:.f90=.o)

all: $(TARGET)

$(TARGET): $(OBJ)
	$(FC) $(LDFLAGS) -o $@ $^

%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@

clean:
	rm -f $(OBJ) $(TARGET) *.mod

data_random:
	$(FC) $(FFLAGS) -o gen_random gen_random.f90
	./gen_random

data_hilbert:
	$(FC) $(FFLAGS) -o gen_hilbert gen_hilbert.f90
	./gen_hilbert

rebuild: clean all

.PHONY: all clean rebuild data_random data_hilbert
