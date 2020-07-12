AC_INIT(mocca, 1.2)

AM_INIT_AUTOMAKE

AC_ARG_ENABLE([shogun],
[  --enable-shogun    Use Shogun machine learning library],
[case "${enableval}" in
  yes) shogun=true ;;
  no)  shogun=false ;;
  *) AC_MSG_ERROR([bad value ${enableval} for --enable-shogun]) ;;
esac],[shogun=false])
AM_CONDITIONAL([USE_SHOGUN], [test x$shogun = xtrue])

AC_PROG_MAKE_SET

AC_LANG(C++)
AC_PROG_CXX

CXXFLAGS="-std=c++11 -fPIC"
LDFLAGS="-pthread"

AC_CONFIG_FILES([Makefile src/Makefile])
AC_OUTPUT

