#nodename=`uname -n`
if [ -e "/usr/bin/g++-4.6" ]; then
    CC=g++-4.6
elif [ -e "/usr/bin/g++-4.5" ]; then
    CC=g++-4.5
fi

CFLAGS="-O3 -march=native -g -std=gnu++0x -Wall -Werror -Wno-deprecated -D_GNU_SOURCE "

echo "$CC $CFLAGS $INC \$@" > $3
chmod +x $3
