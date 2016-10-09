#!/bin/sh
mvn package && mvn install:install-file -Dfile=./target/jbullet-20161008.jar -DgroupId=ca.gsimard.jbullet -DartifactId=jbullet -Dversion=20161008 -Dpackaging=jar -DlocalRepositoryPath=../repo/ && cp -v target/jbullet-20161008.jar ../scripts/server/deps/jbullet.jar && scp target/jbullet-20161008.jar simard@192.168.0.150:./spacecraft/deps/jbullet.jar
