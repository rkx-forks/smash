DEPS="tests.o simplex.o input.o chain.o"
redo-ifchange $DEPS
./compile -L/home/rsbowman/src/boost_1_48_0/stage/lib -o $3 $DEPS -lboost_unit_test_framework
