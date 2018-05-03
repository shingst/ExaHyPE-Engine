#!/bin/bash

#module load intel/16.0
#module load git
#module load cmake
#module load netcdf

local_dir=$HOME/local
temp_build_dir=$HOME/easi_build_dir_mac
mkdir -p ${temp_build_dir}
mkdir -p ${local_dir}

c_compiler=mpicc
cpp_compiler=mpiCC
fortran_compiler=mpif90

cmake_set_install_path="-DCMAKE_INSTALL_PREFIX=${local_dir}"
cmake_set_compiler="-DCMAKE_CXX_COMPILER=${cpp_compiler} -DCMAKE_C_COMPILER=${c_compiler}"

cmake_flags="${cmake_set_install_path} ${cmake_set_compiler}"


export PATH=${local_dir}/bin:$PATH
export CPATH=${local_dir}/include:$CPATH
export CPPPATH=${local_dir}/include:$CPPPATH
export LIBRARY_PATH=${local_dir}/lib:$LIBRARY_PATH
export LD_LIBRARY_PATH=${local_dir}/lib:$LD_LIBRARY_PATH
export PKG_CONFIG_PATH=${local_dir}/lib/pkgconfig:$PKG_CONFIG_PATH


function clone_or_pull {
    if [ ! -d ${1} ]; then
	git clone --recursive ${2} ${1}
    else
	cd ${1}
	git clean -f
	git pull
	cd ..
    fi
}


#check for help2man if not install (required by libtool)
help2man_repo="https://github.com/Distrotech/help2man.git"
help2man_src="$temp_build_dir/help2man"

if [ hash libtool &> /dev/null ]; then
    clone_or_pull ${help2man_src} ${help2man_repo}
    cd ${help2man_src}
    ./configure --prefix=${local_dir}
    make -j8
    make install
fi



#check for libtool if not install (required by current version of numa)
libtool_repo="git://git.savannah.gnu.org/libtool.git"
libtool_src="${temp_build_dir}/libtool"
if [ hash libtool &> /dev/null ]; then
    clone_or_pull ${libtool_src} ${libtool_repo}
    cd ${libtool_src}
    ./bootstrap
    ./configure --prefix=${local_dir}
    make -j8
    make install
fi
export LIBTOOL="libtool"

#clone and build numa
numa_repo="https://github.com/numactl/numactl.git"
numa_build="${temp_build_dir}/numa_build"
numa_src="${temp_build_dir}/numa"

if [ ! -e "${local_dir}/lib/libnuma.so" ]; then
    echo "Building Numa"
    clone_or_pull ${numa_src} ${numa_repo} 
    cd ${numa_src}
    libtoolize
    ./autogen.sh
    ./configure --prefix=${local_dir}
    make -j8
    make install
fi

#clone and build impalaJIT
impala_repo="https://github.com/Manuel1605/ImpalaJIT.git"
impala_build="${temp_build_dir}/ImpalaJIT_build"
impala_src="${temp_build_dir}/ImpalaJIT"

if [ ! -e "${local_dir}/lib/libimpalajit.a" ]; then
    echo "Building Impala"
    clone_or_pull ${impala_src} ${impala_repo}
    mkdir -p ${impala_build}
    cd ${impala_build}
    cmake ${cmake_flags} ${impala_src}
    make -j8
    make install
fi

#clone and build yaml
yaml_repo="https://github.com/jbeder/yaml-cpp.git"
yaml_build="${temp_build_dir}/yaml_build"
yaml_src="${temp_build_dir}/yaml-cpp"

if [ ! -e "${local_dir}/lib/libyaml-cpp.a" ]; then
    echo "Building YAML"
    clone_or_pull ${yaml_src} ${yaml_repo}
    mkdir -p ${yaml_build}
    cd ${yaml_build}
    cmake_flags_yaml="${cmake_flags} -DCMAKE_CXX_FLAGS=--std=c++11"
    cmake ${cmake_flags_yaml} ${yaml_src}
    make -j8
    make install
fi

#build HDF5
hdf5_repo="https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.8/hdf5-1.8.20/src/hdf5-1.8.20.tar"
hdf5_tar=`echo ${hdf5_repo} | awk -F "/" '{print $NF}'`
hdf5_src="${temp_build_dir}/hdf5_src"

if [ ! -e "${local_dir}/lib/libhdf5.a" ]; then
    echo "Install Hdf5"
    mkdir -p ${hdf5_src}
    cd ${hdf5_src}

    if [ ! -e ${hdf5_tar} ]; then
	wget ${hdf5_repo}
    fi
    tmp_dir=`echo ${hdf5_tar} | sed "s/.tar//g" | sed "s/.gz//g" `

    if [ ! -e ${tmp_dir} ]; then
	tar -xvf ${hdf5_tar}
    fi

    cd ${tmp_dir}
    make clean
   CC=${c_compiler} FC=${fortran_compiler} ./configure --enable-parallel --prefix=${local_dir} --with-zlib --disable-shared --enable-fortran
    #knl only !
#    CFLAGS=-fPIC FFLAGS=-fPIC CC=${c_compiler} FC=${fortran_compiler} ./configure --enable-parallel --prefix=${local_dir} --with-zlib --disable-shared --enable-fortran
    make -j8
    make install
fi

netCDF_repo="ftp://ftp.unidata.ucar.edu/pub/netcdf/netcdf-4.5.0.tar.gz"
netCDF_tar=`echo ${netCDF_repo} | awk -F "/" '{print $NF}'`
netCDF_src="${temp_build_dir}/netCDF_src"

if [ ! -e "${local_dir}/lib/libnetcdf.a" ]; then
    echo "Install netCDF"
    mkdir -p ${netCDF_src}
    cd ${netCFD_src}

    if [ ! -e ${netCDF_tar} ]; then
	wget ${netCDF_repo}
    fi

    tmp_dir=`echo ${netCDF_tar} | sed "s/.tar//g" | sed "s/.gz//g" `

    if [ ! -e ${tmp_dir} ]; then
	tar -xvf ${netCDF_tar}
    fi

    cd ${tmp_dir}
    make clean

#    CC=mpiicc CFLAGS=-fPIC ./configure --enable-shared=no --prefix=${local_dir}
    CC=mpiicc ./configure --enable-shared=no --prefix=${local_dir}
    make -j8
    make install
fi

#clone and build ASAGI
asagi_repo="https://github.com/TUM-I5/ASAGI.git"
asagi_build="${temp_build_dir}/asagi_build"
asagi_src="${temp_build_dir}/ASAGI"

if [ ! -e "${local_dir}/lib/libasagi.so" ]; then
    echo "Building ASAGI"
    clone_or_pull ${asagi_src} ${asagi_repo}
    mkdir -p ${asagi_build}
    cd ${asagi_build}
    cmake_flags_asagi="-DNOMPI=OFF -DNONUMA=OFF -DNUMA_LIBRARY=${local_dir}/lib -DNUMA_INCLUDE_DIR=${local_dir}/include ${cmake_flags}"
    cmake ${cmake_flags_asagi}  ${asagi_src}
    make -j8
    make install
    echo "ASAGI Build"
fi

#clone easi
echo "Cloning/Pulling easi"
easi_repo="https://github.com/SeisSol/easi.git"
easi_src="${temp_build_dir}/easi"
easi_build="${temp_build_dir}/easi_build"
clone_or_pull $easi_src $easi_repo
cp -r $easi_src/include/* $local_dir/include/.
rm -rf $easi_build
mkdir $easi_build
cd $easi_build
cmake ${easi_src}

echo "add this to your .bashrc:"
echo 'export PATH=${local_dir}/bin:\$PATH'
echo 'export CPATH=${local_dir}/include:\$CPATH'
echo 'export CPPPATH=${local_dir}/include:\$CPPPATH'
echo 'export LIBRARY_PATH=${local_dir}/lib:\$LIBRARY_PATH'
echo 'export LD_LIBRARY_PATH=${local_dir}/lib:\$LD_LIBRARY_PATH'
echo 'export PKG_CONFIG_PATH=${local_dir}/lib/pkgconfig:\$PKG_CONFIG_PATH'

echo "Done."

cd $HOME

