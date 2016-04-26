#!/bin/bash
retool='samtools-0.1.19'
mkdir -p $PROGRAMS

if [[ $retool = '' ]]; then
	echo 'Do you forgot to add a tool name as parameter? Or'
	read -p 'do you want to recompile all tools [y|n] ' in
	if [[ $in = 'y' ]] || [[ $in = 'Y' ]]; then		
		retool='all'
	else 
		read -p 'Please enter the correct tool name as used in '$PROGRAMS'/<toolname-version> or hit return to exit' retool
		if [[ "$retool" = '' ]]; then
			echo 'Recompilation aborted...'
			exit
		fi
	fi
fi

pwd=$PWD
bit=$(uname -m | awk '{if ($0 == "i686"){print "32-"}else{print ""}}')
if [[ $OSTYPE == darwin* ]]; then 
	bit='mac-'
fi

download () {
	if [[ ! -d $PROGRAMS/$tool ]] && [[ $ex -eq 0 ]]; then
		echo 'Downloading '$tool	
		wget -T 10 -N www.rna.uni-jena.de/supplements/gorap/$bit$tool.tar.gz
		if [[ $? -gt 0 ]]; then
			echo $tool' download failed'
			exit 1
		fi 
		echo 'Extracting '$tool
		tar -xzf $bit$tool.tar.gz -C $PROGRAMS
		rm $bit$tool.tar.gz
	fi
}

tool='samtools-0.1.19'
samtool=$tool
ex=0
download

tool='zlib-1.2.8'		
ex=0
download				

cd $PROGRAMS/$tool 
make clean
./configure --prefix=$PROGRAMS/$samtool				
make
if [[ $? -gt 0 ]]; then				
	exit 1
fi
make install
make clean

tool='ncurses-5.9'
ex=0
download

cd $PROGRAMS/$tool  
make clean
./configure --prefix=$PROGRAMS/$samtool
make
if [[ $? -gt 0 ]]; then				
	exit 1
fi
make install
make clean

cd $PROGRAMS/$samtool
sed -i '/^CFLAGS/ c\CFLAGS = -g -Wall -O2 -fPIC -I. -I.. -I./include -I../include -L. -L.. -L./lib -L../lib #-m64 #-arch ppc' Makefile
make clean
cp -r include/* .
cp -r lib/* .                
make 
if [[ $? -gt 0 ]]; then				
	exit 1
fi                       
mkdir -p bin 
cp samtools bin/samtools
cd $pwd
