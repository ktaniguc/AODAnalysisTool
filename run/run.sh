#!/bin/bash
#OUTPUT_DIR="$TestArea/output"
USE_BSUB="1" ### 0: no use lsf batch  1: use lsf batch
SEPARATE="1" ### 0: use all files in 1 job 1: one job per one file( only work with batch job )
#BSUB_QUE="30m"
#BSUB_QUE="12h"
BSUB_QUE="1d"
GRL_NAME="current_grl.xml"
DATE=$(date '+%Y_%m_%d_%H_%M')

if [ "$TestArea" = "" ] ; then
  echo "TestArea is missing"
  exit 1
fi

source $TestArea/x86_64-centos7-gcc8-opt/setup.sh

if [ $USE_BSUB = 1 ] ; then
  echo "Batch will be used..." && sleep 1
else
  echo "local job will be used..." && sleep 2
fi

if [ "$1" = "" ] ; then
  if [ "$USE_BSUB" = "1" ] ; then
    bsub -q $BSUB_QUE -o ./LOGS/log_${DATE}.txt athena.py CalcEffAlg_optionsMT.py
  else
    athena.py CalcEffAlg_optionsMT.py &> log_1101_01 &
  fi
  exit 0
fi

LIST_NAME=$1
LOG_INDEX="test"
LOG_NAME="log_"$LOG_INDEX"_"$DATE
number=1

if [ -f "$LIST_NAME" ] ; then
  echo $LIST_NAME
  if [ "$USE_BSUB" = "1" ] ; then
    if [ "$SEPARATE" = "0" ] ; then
      echo "athena.py -c 'LocalInputFileList=\"$LIST_NAME\"' CalcEffAlg_optionsMT.py" > tmp_athena_run.sh
      cat tmp_athena_run.sh
      bsub -q $BSUB_QUE -o log_1101 < tmp_athena_run.sh
    else
      for FILE in `cat $LIST_NAME`
      do
        echo ""
        DIR="Output/$DATE"
        DIR_TMP="Output/$DATE/FILE_${number}"
        mkdir -p $DIR_TMP
        echo "$FILE" > $DIR_TMP/tmp.list
        cp "$GRL_NAME" $DIR_TMP/
        cp $LIST_NAME $DIR_TMP/
        cp CalcEffAlg_optionsMT.py $DIR_TMP/
        cp run.sh $DIR_TMP/
        cd $DIR_TMP
        TMP_LIST="$DIR_TMP/tmp.list"
        echo "athena.py -c 'LocalInputFileList=\"tmp.list\"' CalcEffAlg_optionsMT.py" > tmp_athena_run.sh
        echo "DIR: "$DIR_TMP
        cat tmp_athena_run.sh
        bsub -q $BSUB_QUE -o $LOG_NAME < tmp_athena_run.sh
        number=`expr $number + 1`
        cd -
      done
    fi
  else
    echo "athena.py -c 'LocalInputFileList=\"$LIST_NAME\"' CalcEffAlg_optionsMT.py &> log1 &" > tmp_athena_run.sh
    cat tmp_athena_run.sh
    sh tmp_athena_run.sh
  fi
  #athena.py -c "LocalInputFileList=\"$LIST_NAME\"" CalcEffAlg_optionsMT.py &> log &
fi
