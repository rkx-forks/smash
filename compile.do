#nodename=`uname -n`
if [ -e "/usr/bin/g++-4.9" ]; then
    CC=g++-4.9
    STD=c++14
elif [ -e "/usr/bin/g++-4.6" ]; then
    CC=g++-4.6
    STD=c++0x
fi
DEBUG_FLAGS="-g"
LIBS="-lboost_program_options -lboost_timer"
CFLAGS="-O1 -march=native --std=$STD -Wall -Werror -Wno-deprecated -D_GNU_SOURCE $LIBS"

echo "$CC $DEBUG_FLAGS $CFLAGS $INC \$@" > $3
chmod +x $3
