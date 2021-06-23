

all: MOCCA


MOCCA:
	emcc ../src/vaux.cpp ../src/main.cpp ../src/sequencelist.cpp ../src/validation.cpp ../src/config.cpp ../src/motifs.cpp ../src/sequences.cpp ../src/models/baseclassifier.cpp ../src/models/seqlo.cpp ../src/models/cpredictor.cpp ../src/models/seqrf.cpp ../src/models/dummypredictor.cpp ../src/models/seqsvm.cpp ../src/models/features.cpp ../src/models/sequenceclassifier.cpp ../src/models/rfmocca.cpp ../src/models/svmmocca.cpp ../src/models/seqdummy.cpp ../src/lib/libsvm-3.17/svm.cpp ../src/lib/ranger/src/Data.cpp ../src/lib/ranger/src/Forest.cpp ../src/lib/ranger/src/ForestClassification.cpp ../src/lib/ranger/src/ForestProbability.cpp ../src/lib/ranger/src/ForestRegression.cpp ../src/lib/ranger/src/Tree.cpp ../src/lib/ranger/src/TreeClassification.cpp ../src/lib/ranger/src/TreeProbability.cpp ../src/lib/ranger/src/TreeRegression.cpp ../src/lib/ranger/src/utility.cpp -o mocca.html -O3 -s EXPORTED_FUNCTIONS='["_main"]' -s FORCE_FILESYSTEM=1 -s ALLOW_MEMORY_GROWTH=1 -s ASYNCIFY -s EXTRA_EXPORTED_RUNTIME_METHODS=FS -s NO_EXIT_RUNTIME=1

serve:
	python3 -m http.server 80


