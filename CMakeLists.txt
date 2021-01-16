cmake_minimum_required(VERSION 3.5)

project (mocca)

set(CMAKE_CXX_STANDARD 11)
set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)

add_library(libsvm STATIC 
    src/lib/libsvm-3.17/svm.cpp
)

add_library(ranger STATIC 
    src/lib/ranger/src/Data.cpp
    src/lib/ranger/src/Forest.cpp
    src/lib/ranger/src/ForestClassification.cpp
    src/lib/ranger/src/ForestProbability.cpp
    src/lib/ranger/src/ForestRegression.cpp
    src/lib/ranger/src/Tree.cpp
    src/lib/ranger/src/TreeClassification.cpp
    src/lib/ranger/src/TreeProbability.cpp
    src/lib/ranger/src/TreeRegression.cpp
    src/lib/ranger/src/utility.cpp
)

set(SOURCES
    src/main.cpp
    src/aux.cpp
    src/config.cpp
    src/motifs.cpp
    src/sequencelist.cpp
    src/sequences.cpp
    src/validation.cpp
    src/models/features.cpp
    src/models/sequenceclassifier.cpp
    src/models/baseclassifier.cpp
    src/models/seqdummy.cpp
    src/models/seqlo.cpp
    src/models/seqsvm.cpp
    src/models/seqrf.cpp
    src/models/dummypredictor.cpp
    src/models/cpredictor.cpp
    src/models/svmmocca.cpp
    src/models/rfmocca.cpp
)

set(LIBS
    libsvm
    ranger
    Threads::Threads
)

OPTION(USE_SHOGUN "Use Shogun." OFF)
IF(USE_SHOGUN)
	set(SOURCES
		${SOURCES}
		src/models/seqlda.cpp
		src/models/seqperceptron.cpp
	)
	set(LIBS
		${LIBS}
		shogun
	)
	add_definitions(-DUSE_SHOGUN)
ENDIF(USE_SHOGUN)

add_executable(mocca ${SOURCES})

target_link_libraries(mocca PRIVATE ${LIBS})

install (TARGETS mocca DESTINATION bin)

#set(CPACK_GENERATOR "DEB")
#set(CPACK_PACKAGE_DESCRIPTION "MOCCA suite")
#set(CPACK_PACKAGE_CONTACT "Bjørn Bredesen bjorn@bjornbredesen.no")
#set(CPACK_PACKAGE_VERSION_MAJOR "1")
#set(CPACK_PACKAGE_VERSION_MINOR "3")
#set(CPACK_PACKAGE_VERSION_PATCH "0")
#include(CPack)

