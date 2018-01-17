TARGET	= Part
OBJS_A	= $(TARGET)_A.o
OBJS_B  = $(TARGET)_B.o
OBJS_C  = $(TARGET)_C.o
TARGET_ext = .prog

FC  	= mpif90
FFLAGS 	= -free

all: $(TARGET)

$(TARGET)_A: $(OBJS_A)
	$(FC) $(FFLAGS) $(OBJS_A) -o $(TARGET)_A$(TARGET_ext)

$(TARGET)_B: $(OBJS_B)
	$(FC) $(FFLAGS) $(OBJS_B) -o $(TARGET)_B$(TARGET_ext)

$(TARGET)_C: $(OBJS_C)
	$(FC) $(FFLAGS) $(OBJS_C) -o $(TARGET)_C$(TARGET_ext)

clean:
	@/bin/rm -fv $(TARGET)_A$(TARGET_ext)
	@/bin/rm -fv $(TARGET)_B$(TARGET_ext)
	@/bin/rm -fv $(TARGET)_C$(TARGET_ext)
	@/bin/rm -fv $(OBJS_A)
	@/bin/rm -fv $(OBJS_B)
	@/bin/rm -fv $(OBJS_C)

%.o: %.f90
	$(FC) $(FFLAGS) -c $<

# dependencies
$(TARGET)_A.o : $(TARGET)_A.f90
$(TARGET)_B.o : $(TARGET)_B.f90
$(TARGET)_C.o : $(TARGET)_C.f90


