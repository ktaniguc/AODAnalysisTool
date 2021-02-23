
echo ""
echo "TestArea = "$TestArea
echo ""


#dir=$(cd $(dirname ${BASH_SOURCE:-$0}); pwd)
dir=$PWD
version=( `echo $AtlasVersion | tr -s '.' ' '` )
project=$AtlasProject

if [ ${version[0]} -le 20 ]; then
    echo "Release" ${version[0]} "=> CMT build"
    cd ${TestArea}
    cmt find_packages
    cmt compile
    cd -


else
    echo "Release" ${version[0]} "=> CMake build"
  
    cd $dir/build/
    if [ $# -gt 0 ] && [ $1 = "cmake" ]; then
        if [ `echo "$AtlasProject" | grep "Athena"` ];then
            echo "cmake " $dir/source/Projects/WorkDir
            cmake -DATLAS_PACKAGE_FILTER_FILE=$dir/package_filters.txt $dir/source/Projects/WorkDir
        else
            echo "cmake " $dir/source/
            cmake $dir/source/
        fi
    fi
    if [ $# -gt 0 ] && [ $1 = "clean" ]; then
        if [ `echo "$AtlasProject" | grep "Athena"` ];then
	    rm -r $dir/build/*
            echo "cmake " $dir/source/Projects/WorkDir
            cmake -DATLAS_PACKAGE_FILTER_FILE=$dir/package_filters.txt $dir/source/Projects/WorkDir
        else
	    rm -r $dir/build/*
            echo "cmake " $dir/source/
            cmake $dir/source/
        fi
    fi
    make -j8
    cd -
fi

#cd $TestArea/WorkArea/cmt
#cmt br cmt config
#cmt br cmt make
#cd -


echo ""
