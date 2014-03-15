#! /bin/sh

SCRIPT_DIR=`dirname $0`
PEPR_HOME=$SCRIPT_DIR/..
LIB_DIR=$PEPR_HOME/lib/
CLASSPATH=.
CLASSPATH=$CLASSPATH:$PEPR_HOME/lib/pepr.jar
CLASSPATH=$CLASSPATH:$PEPR_HOME/lib/log4j.jar

echo $CLASSPATH

MAINCLASS=edu.vt.vbi.ci.pepr.tree.pipeline.PhyloPipeline 

set -x
JAVA_CMD="java -Xmx4g -Dlog4j.configuration=file:$PEPR_HOME/lib/log4j.properties\
           -cp $CLASSPATH $MAINCLASS $*"

$JAVA_CMD

