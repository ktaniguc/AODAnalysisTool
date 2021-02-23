
#dir=$(cd $(dirname ${BASH_SOURCE:-$0}); pwd)
setupATLAS

lsetup git

BUILD_DIR=$PWD/build
#rm -r $BUILD_DIR
if [ -d "$BUILD_DIR" ] ; then
  echo "Directory for build is $BUILD_DIR"
else
  echo "Creating $BUILD_DIR ..." && mkdir -p $BUILD_DIR
fi
cd $BUILD_DIR
#asetup Athena,22.0.0
#asetup Athena,master,latest,here 
asetup Athena,22.0.24,here
#asetup Athena,22.0.15,here
#asetup Athena,22.0.8,here
#asetup Athena,21.0.90,here
cd -
