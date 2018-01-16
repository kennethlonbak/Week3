TARGET	= Hallo_World_MPI
OBJS	= $(TARGET).o
TARGET_ext = .prog

FC  	= mpif90
#FFLAGS 	= -free -ftrap=%none
#FC 		= gfortran
FFLAGS 	= -free
all: $(TARGET)

$(TARGET): $(OBJS)
	$(FC) $(FFLAGS) $(OBJS) -o $(TARGET)$(TARGET_ext)

clean:
	@/bin/rm -fv $(TARGET)$(TARGET_ext)
	@/bin/rm -fv $(OBJS)


%.o: %.f90
	$(FC) $(FFLAGS) -c $<

# dependencies
$(TARGET).o : $(TARGET).f90


