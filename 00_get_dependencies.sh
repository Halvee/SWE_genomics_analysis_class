#!/bin/bash

# specific to mac users :
# 1. need to install brew
# (https://www.datacamp.com/tutorial/homebrew-install-use)
# /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
# 2. need to install git
# brew install git

# make directories for ..
# 1. where the binaries will go
mkdir -p bin/
# 2. where some of the non-binary source code will go
mkdir -p src/
# 3. where the data will go
mkdir -p data/

# get binaries (mainly PLINK)
# 1. go into binary dir
cd bin/
# 2. download plink 1.9 (note for mac version if you're using M1 then you need Rosetta 2 installed)
wget -N "https://s3.amazonaws.com/plink1-assets/plink_mac_20230116.zip"
# wget -N "https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20230116.zip"
# wget -N "https://s3.amazonaws.com/plink1-assets/plink_win64_20230116.zip"
# 3. download plink 2.0
wget -N "https://s3.amazonaws.com/plink2-assets/plink2_mac_arm64_20230915.zip"
# wget -N "https://s3.amazonaws.com/plink2-assets/plink2_mac_20230915.zip"
# wget -N "https://s3.amazonaws.com/plink2-assets/plink2_linux_x86_64_20230915.zip"
# 4. unzip plink 1.9 and plink 2.0
unzip plink_mac_20230116.zip; unzip plink2_mac_arm64_20230915.zip
# unzip plink_linux_x86_64_20230116.zip; unzip plink2_linux_x86_64_20230915.zip
# unzip plink_mac_20230116.zip; unzip plink2_mac_20230915.zip
# 5. install GCTA
wget -N "https://yanglab.westlake.edu.cn/software/gcta/bin/gcta-1.94.2-MacOS-ARM-x86_64.zip"
# wget -N "https://yanglab.westlake.edu.cn/software/gcta/bin/gcta-1.94.1-Win-x86_64.zip"
# wget -N "https://yanglab.westlake.edu.cn/software/gcta/bin/gcta-1.94.1-linux-kernel-3-x86_64.zip"
unzip gcta-1.94.2-MacOS-ARM-x86_64.zip
# unzip gcta-1.94.1-Win-x86_64.zip
# unzip gcta-1.94.1-linux-kernel-3-x86_64.zip
# go back to parent directory
cd ../

# lets download the data we'll be using. 
# 1. go into the data directory
cd data/
# 2. get the 'download all' master list for synthetic gwas data using wget
wget -N "https://personal.broadinstitute.org/sawasthi/share_links/UzoZK7Yfd7nTzIxHamCh1rSOiIOSdj_gwas-qcerrors.py/download_all.txt"
# 3. make list of properly formatted file URLs and then download with wget
cut -d" " -f 2 download_all.txt > files_download.list
wget -N -i files_download.list
# 4. extra credit : let's check to see if the file checksums equal those reported (for data integrity)
#    and put the results of the checksum check into a single file.
> file_cksum.check.txt
for file in *.bed *.bim *.fam
do 
  cksum1=`cksum $file`
  cksum2=`cat $file.cksum.txt`
  if [[ $cksum1 == $cksum2 ]]
  then
    echo "$file : PASS" >> file_cksum.check.txt
  else
    echo "$file : FAIL" >> file_cksum.check.txt
  fi
done
# download PLINK formatted HGDP files
wget -N "https://www.dropbox.com/s/hppj1g1gzygcocq/hgdp_all.pgen.zst" \
        "https://www.dropbox.com/s/1mmkq0bd9ax8rng/hgdp_all.pvar.zst" \
        "https://www.dropbox.com/s/0zg57558fqpj3w1/hgdp.psam"
# go back to parent directory
cd ../

# lets install a few programs with source code into the src dir
# 1. go into src dir
cd src/
# 2. install ldsc
git clone https://github.com/bulik/ldsc.git 
# 3. install PRS-CS
git clone https://github.com/getian107/PRScs.git
# 4. install peddy
git clone https://github.com/brentp/peddy.git
# go back to parent dir
cd ../

# other resources : 
# datacamp.com - free courses for learning R, python, etc
# https://www.freecodecamp.org/ - free courses for R, python with more direct links to youtube vids of walkthroughs
# coursera
#   hands-on intro to linux commands and shell scripting:
#   https://www.coursera.org/learn/hands-on-introduction-to-linux-commands-and-shell-scripting
