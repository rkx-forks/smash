DEPS="main_collapse.o simplex.o input.o chain.o"
redo-ifchange compile $DEPS
./compile -o $3 $DEPS -lboost_program_options -lboost_timer -lboost_system
