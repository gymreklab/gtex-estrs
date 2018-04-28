AWS_ACCESS_KEY=$1
AWS_SECRET_KEY=$2

die()
{
    BASE=$(basename -- "$0")
    echo "$BASE: error: $@" >&2
    terminate
    exit 1
}

HOMEDIR=/root/

AWS_DIR=${HOMEDIR}/.aws
AWS_CONFIG_FILE=${AWS_DIR}/config
AWS_CRED_FILE=${AWS_DIR}/credentials

# Install things
sudo apt-get update || die "Could not update"
sudo apt-get -y install awscli || die "Could not install aws"
sudo apt-get -y install git || die "Could not install git"
sudo apt-get -y install make gcc libz-dev libncurses5-dev libbz2-dev liblzma-dev libcurl3-dev libssl-dev autoconf || die "Could not install devtools"
sudo apt-get install -y python-setuptools python-dev build-essential || die "Could not install python"
sudo easy_install pip || die "Could not install pip"
sudo pip install joblib numpy scipy pytabix pyvcf statsmodels pycallgraph pysam || die "Could not install python libraries"

mkdir -p ${HOMEDIR}/source || die "Could not create source dir"
cd ${HOMEDIR}/source || die "Could not go to source dir"

cd ${HOMEDIR}/source
wget https://github.com/samtools/bcftools/releases/download/1.4.1/bcftools-1.4.1.tar.bz2
tar -xvf bcftools-1.4.1.tar.bz2
cd bcftools-1.4.1
make
sudo make install

cd ${HOMEDIR}
git clone https://github.com/samtools/htslib
cd htslib
autoheader
autoconf
./configure --enable-libcurl
make
sudo make install

cd ${HOMEDIR}
git clone https://github.com/samtools/samtools
cd samtools
autoconf -Wno-syntax 
./configure
make
sudo make install

cd ${HOMEDIR}
git clone https://github.com/tfwillems/HipSTR
cd HipSTR
make
cp HipSTR /usr/local/bin

# Set up AWS credentials
echo "Setting up AWS credentials in ${AWS_DIR}"
mkdir -p ${AWS_DIR} || die "Could not create AWS dir"
echo "[default]" > ${AWS_CONFIG_FILE} || die "Could not write to ${AWS_CONFIG_FILE}"
echo "output = table" >> ${AWS_CONFIG_FILE} || die "Could not write to ${AWS_CONFIG_FILE}"
echo "region = us-east-1" >> ${AWS_CONFIG_FILE}  || die "Could not write to ${AWS_CONFIG_FILE}"
echo "aws_secret_access_key = ${AWS_SECRET_KEY}" >> ${AWS_CONFIG_FILE} || die "Could not write to ${AWS_CRED_FILE}"
echo "aws_access_key_id = ${AWS_ACCESS_KEY}" >> ${AWS_CONFIG_FILE} || die "Could not write to ${AWS_CRED_FILE}"
echo "[default]" > ${AWS_CRED_FILE} || die "Could not write to ${AWS_CRED_FILE}"
echo "aws_secret_access_key = ${AWS_SECRET_KEY}" >> ${AWS_CRED_FILE} || die "Could not write to ${AWS_CRED_FILE}"
echo "aws_access_key_id = ${AWS_ACCESS_KEY}" >> ${AWS_CRED_FILE} || die "Could not write to ${AWS_CRED_FILE}"

# SRA toolkit
cd ${HOMEDIR}/source
wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.8.2-1/sratoolkit.2.8.2-1-ubuntu64.tar.gz
tar -xzvf sratoolkit.2.8.2-1-ubuntu64.tar.gz
cp sratoolkit.2.8.2-1-ubuntu64/bin/* /usr/local/bin/

# Set up dbgap credentials in /root/dbgap (manusally)
# Set up using vdb-config (manually)
sudo mkdir /storage/gtex-data

# Install apsera
cd ${HOMEDIR}/source
wget http://download.asperasoft.com/download/sw/connect/3.7.4/aspera-connect-3.7.4.147727-linux-64.tar.gz
tar -xzvf aspera-connect-3.7.4.147727-linux-64.tar.gz 
# Manually install as ubuntu user, can'tinstall as root
sudo cp /home/ubuntu/.aspera/connect/bin/ascp /usr/local/bin/

# Download what need for hipstr
#scp -i ~/keys/micro_key.pem /storage/resources/dbase/human/hg19/Homo_sapiens_assembly19.fasta ubuntu@ec2-34-224-22-182.compute-1.amazonaws.com:/mnt/resources/
sudo mkdir /mnt/resources
sudo chmod 777 /mnt/resources
cd /mnt/resources/
wget https://github.com/HipSTR-Tool/HipSTR-references/raw/master/human/GRCh37.hipstr_reference.bed.gz
gunzip GRCh37.hipstr_reference.bed.gz

# Get github for gtex project
cd ${HOMEDIR}/source
git clone https://github.com/gymreklab/gtex-estrs
