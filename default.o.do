redo-ifchange compile $2.cpp
./compile -MD -MF $2.d -c -o $3 $2.cpp
read DEPS <$2.d
redo-ifchange ${DEPS#*:}
rm $2.d
