#! /bin/sh

PEPRLIB=${0/scripts\/pepr.sh/lib/}
PEPRLIB=${PEPRLIB/.\/pepr.sh/..\/lib/}
PEPRLIB=${PEPRLIB/pepr.sh/..\/lib/}

CLASSPATH=.
CLASSPATH=$CLASSPATH:$PEPRLIB/pepr.jar
CLASSPATH=$CLASSPATH:$PEPRLIB/log4j.jar

MAINCLASS=edu.vt.vbi.ci.pepr.tree.pipeline.PhyloPipeline

set -x

java -Xmx8g -Djava.util.Arrays.useLegacyMergeSort=true -Dlogfile.name=pepr.log -Dlog4j.configuration=file:$PEPRLIB/log4j.properties -cp $CLASSPATH $MAINCLASS $*
