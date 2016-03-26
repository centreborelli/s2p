#!/bin/sh -e
#
# Check there is no memory leak with valgrind/memcheck.

_test_memcheck() {
    test "0" = "$( valgrind --tool=memcheck $* 2>&1 | grep -c 'LEAK' )"
}

################################################

_log_init

echo "* check memory leaks"
_log make distclean
_log make
_log _test_memcheck example/rand
_log _test_memcheck example/rand 0
_log make distclean
_log make CFLAGS=
_log _test_memcheck example/rand
_log _test_memcheck example/rand 0

_log make distclean

_log_clean
