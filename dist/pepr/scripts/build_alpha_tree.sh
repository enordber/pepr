#! /bin/sh

PEPRLIB=${0/scripts\/pepr.sh/lib/}
PEPRLIB=${PEPRLIB/.\/pepr.sh/..\/lib/}
PEPRLIB=${PEPRLIB/pepr.sh/..\/lib/}

ALPHADIR=$PEPRLIB"../Alphaproteobacteria/"
GENOMEFILES=$ALPHADIR"*.faa"
OUTGROUPFILES=$ALPHADIR"/outgroup/*faa"
CLASSPATH=.
CLASSPATH=$CLASSPATH:$PEPRLIB/pepr.jar
CLASSPATH=$CLASSPATH:$PEPRLIB/log4j.jar

OPTIONS=" -track blast_raxml"
OPTIONS=$OPTIONS" -genome "$GENOMEFILES
OPTIONS=$OPTIONS" -outgroup "$OUTGROUPFILES
OPTIONS=$OPTIONS" -outgroup_count 8"
OPTIONS=$OPTIONS" -ml_matrix PROTGAMMAWAG"

MAINCLASS=edu.vt.vbi.ci.pepr.tree.pipeline.PhyloPipeline

set -x

java -Dlog4j.configuration=file:$PEPRLIB/log4j.properties -cp $CLASSPATH $MAINCLASS $OPTIONS $*
