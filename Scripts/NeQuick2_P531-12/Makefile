ALLPROGS = slQu_2 NeQVal
COMMAND = g77
# or f77, g77 or g77_win32 depending on the compiler and platform
include rules_$(COMMAND).make

# or gcc or g77
all: $(ALLPROGS)

slQu_2: \
	NeQuick_2.o \
	slQu_2.for
	$(CC) $(FLAGS) $^ -o $@ $(LIBS) 

NeQVal: \
	NeQuick_2.o \
	NeQVal.for
	$(CC) $(FLAGS) $^ -o $@ $(LIBS)

NeQuick_2.o: \
	NeQuick_2.for
	$(CC) $(XFLAGS) -c $^ -o $@
clean: 
	rm -f *.o $(ALLPROGS) core *~ *.exe
