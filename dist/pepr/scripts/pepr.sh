#! /bin/sh

PEPRLIB=${0/scripts\/pepr.sh/lib/}
PEPRLIB=${PEPRLIB/.\/pepr.sh/..\/lib/}
PEPRLIB=${PEPRLIB/pepr.sh/..\/lib/}

CLASSPATH=.
CLASSPATH=$CLASSPATH:$PEPRLIB/pepr.jar
CLASSPATH=$CLASSPATH:$PEPRLIB/log4j.jar

OPTIONS=" -track blast_raxml"
OPTIONS=$OPTIONS" -outgroup_count 1"
OPTIONS=$OPTIONS" -ml_matrix PROTGAMMAWAG"
OPTIONS=$OPTIONS" -refine true"

MAINCLASS=edu.vt.vbi.ci.pepr.tree.pipeline.PhyloPipeline

set -x

java -Xmx8g -Djava.util.Arrays.useLegacyMergeSort=true -XX:ParallelGCThreads=1 -Dlogfile.name=pepr.log -Dlog4j.configuration=file:$PEPRLIB/log4j.properties -cp $CLASSPATH $MAINCLASS $OPTIONS $*
