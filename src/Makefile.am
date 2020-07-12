AUTOMAKE_OPTIONS = foreign subdir-objects
bin_PROGRAMS = mocca

if USE_SHOGUN
mocca_SOURCES = main.cpp aux.cpp config.cpp motifs.cpp sequencelist.cpp sequences.cpp validation.cpp models/baseclassifier.cpp models/cpredictor.cpp models/dummypredictor.cpp models/features.cpp models/seqdummy.cpp models/seqlo.cpp models/seqsvm.cpp models/seqrf.cpp models/seqlda.cpp models/seqperceptron.cpp models/sequenceclassifier.cpp models/svmmocca.cpp models/rfmocca.cpp lib/libsvm-3.17/svm.cpp lib/ranger/src/Forest.cpp lib/ranger/src/ForestClassification.cpp lib/ranger/src/ForestProbability.cpp lib/ranger/src/ForestRegression.cpp lib/ranger/src/Tree.cpp lib/ranger/src/TreeClassification.cpp lib/ranger/src/TreeProbability.cpp lib/ranger/src/TreeRegression.cpp lib/ranger/src/Data.cpp lib/ranger/src/utility.cpp
mocca_LDFLAGS = -O3
mocca_LDADD = -lshogun
mocca_CPPFLAGS = -DUSE_SHOGUN
else
mocca_SOURCES = main.cpp aux.cpp config.cpp motifs.cpp sequencelist.cpp sequences.cpp validation.cpp models/baseclassifier.cpp models/cpredictor.cpp models/dummypredictor.cpp models/features.cpp models/seqdummy.cpp models/seqlo.cpp models/seqsvm.cpp models/seqrf.cpp models/sequenceclassifier.cpp models/svmmocca.cpp models/rfmocca.cpp lib/libsvm-3.17/svm.cpp lib/ranger/src/Forest.cpp lib/ranger/src/ForestClassification.cpp lib/ranger/src/ForestProbability.cpp lib/ranger/src/ForestRegression.cpp lib/ranger/src/Tree.cpp lib/ranger/src/TreeClassification.cpp lib/ranger/src/TreeProbability.cpp lib/ranger/src/TreeRegression.cpp lib/ranger/src/Data.cpp lib/ranger/src/utility.cpp
mocca_LDFLAGS = -O3
endif

