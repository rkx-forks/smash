#nodename=`uname -n`
if [ -e "/usr/bin/g++-4.9" ]; then
    CC=g++-4.9
    STD=c++14
elif [ -e "/usr/bin/g++-4.6" ]; then
    CC=g++-4.6
    STD=c++1y
fi
CFLAGS="-O3 -march=native --std=$STD -Wall -Werror -Wno-deprecated -D_GNU_SOURCE -lboost_program_options"

echo "$CC $CFLAGS $INC \$@" > $3
chmod +x $3
