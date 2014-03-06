#! /bin/sh

CLASSPATH=.
CLASSPATH=$CLASSPATH:../lib/pepr-20140221.jar
CLASSPATH=$CLASSPATH:../lib/log4j.jar

echo $CLASSPATH

MAINCLASS=edu.vt.vbi.ci.pepr.tree.pipeline.PhyloPipeline 

set +x

java -cp $CLASSPATH $MAINCLASS 
