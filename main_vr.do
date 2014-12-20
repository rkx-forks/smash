DEPS="main_vr.o simplex.o input.o chain.o"
redo-ifchange $DEPS
./compile $CFLAGS $LIBS -o $3 $DEPS
